"""The custom-database build use-case, shared by the CLI and MCP.

Turns a user-supplied input into a fully built gapit-native database: FASTA
inputs run the gene pipeline (records.jsonl -> sequences + BLAST index +
manifest, written last), GBK/GFF inputs run the cluster pipeline (locus
FASTA + features.json + manifest ``kind: cluster``, gapit.dbbuild). This
module is orchestration plus kind dispatch + a metadata merge only — the
building blocks live in records/dbbuild/fasta/dbcodec/db (SPEC.md §11). The
typer command (:mod:`gapit.cmd_db_build`) and the MCP tool
(:mod:`gapit.mcp_tools`) are thin callers.

Header kind is detected PER RECORD (mixed files allowed): a ``gapit|``
prefix decodes through the strict tagged codec, anything else through the
abricate ``~~~`` rules (a plain id carries no ``~~~`` and decodes to itself).
``db`` is always the TARGET name and ``source_id`` keeps the original id
token; record order is input order. Sequences are stored verbatim — no
provider-style normalization, this is the user's curated truth.
"""

import csv
import shutil
from collections.abc import Callable
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

from pydantic import BaseModel

from gapit import config
from gapit.clusterbuild import is_cluster_input, perform_cluster_build
from gapit.db import mol_type
from gapit.dbbuild import build_database
from gapit.dbcodec import decode_seqid
from gapit.errors import DatabaseError, InputError, UsageError
from gapit.fasta import FastaRecord, iter_fasta
from gapit.records import Record, write_records
from gapit.typing_gene import validate_gene_typing
from gapit.typing_models import read_typing_document, typing_schema_of

Dbtype = Literal["nucl", "prot"]
Kind = Literal["gene", "cluster"]


class BuildReceipt(BaseModel, frozen=True):
    """One-line JSON success receipt for `db build` (mirrors db fetch);
    ``records`` counts genes for a gene build, loci for a cluster build."""

    db: str
    records: int
    dbtype: Dbtype
    destination: str


@dataclass(frozen=True, slots=True)
class Metadata:
    """Parsed --tsv: which merge columns the header carries, and the row
    values per gene (accession, ';'-joined function; '' for absent columns)."""

    columns: frozenset[str]
    rows: dict[str, tuple[str, str]]


def _split_function(joined: str) -> tuple[str, ...]:
    """';'-joined classes -> tuple, empty pieces dropped."""
    return tuple(part for part in joined.split(";") if part)


def _to_record(fasta: FastaRecord, name: str, default_product: str) -> Record:
    """One FASTA record -> one Record.

    decode_seqid routes by header kind; the product prefers the FASTA
    description, then --description, then the gene (the provider-side
    ``description or gene`` pattern). Malformed ``gapit|`` headers raise
    their native HEADER_MALFORMED DatabaseError.
    """
    header = decode_seqid(fasta.id, default_db=name)
    return Record(
        db=name,
        gene=header.gene,
        sequence=fasta.sequence,
        accession=header.accession,
        function=_split_function(header.function),
        product=fasta.description or default_product or header.gene,
        source_id=fasta.id,
    )


def _read_metadata(path: Path, warn: Callable[[str], None]) -> Metadata:
    """Parse the metadata TSV: header row mandatory and must carry ``gene``
    (InputError METADATA_MALFORMED otherwise); ``accession``/``function``
    are optional per-file, extra columns are ignored. Duplicate genes keep
    the FIRST row (warn); short rows read their missing cells as empty.
    """
    if not path.is_file():
        raise InputError(
            f"metadata file not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        )
    with path.open(encoding="utf-8", newline="") as handle:
        rows = csv.reader(handle, delimiter="\t")
        header = next((row for row in rows if row), None)
        if header is None or "gene" not in header:
            raise InputError(
                f"metadata TSV needs a header row with a 'gene' column: {path}",
                code="METADATA_MALFORMED",
                context={"file": str(path), "needed": "gene"},
            )
        columns = {column: index for index, column in enumerate(header)}

        def cell(padded: list[str], column: str) -> str:
            return padded[columns[column]]

        metadata: dict[str, tuple[str, str]] = {}
        for row in rows:
            if not any(value.strip() for value in row):
                continue
            padded = row + [""] * (len(header) - len(row))
            gene = cell(padded, "gene")
            if not gene:
                continue
            if gene in metadata:
                warn(f"duplicate gene {gene!r} in metadata TSV: keeping the first row")
                continue
            metadata[gene] = (
                cell(padded, "accession") if "accession" in columns else "",
                cell(padded, "function") if "function" in columns else "",
            )
    return Metadata(columns=frozenset(columns), rows=metadata)


def _merge(records: list[Record], metadata: Metadata, warn: Callable[[str], None]) -> list[Record]:
    """Overwrite accession/function on records whose gene has a metadata row
    (only the columns the TSV actually carries); genes only present in the
    TSV are warned about and skipped."""
    merged: list[Record] = []
    for record in records:
        row = metadata.rows.get(record.gene)
        if row is None:
            merged.append(record)
            continue
        update: dict[str, object] = {}
        if "accession" in metadata.columns:
            update["accession"] = row[0]
        if "function" in metadata.columns:
            update["function"] = _split_function(row[1])
        merged.append(record.model_copy(update=update))
    fasta_genes = {record.gene for record in records}
    for gene in sorted(set(metadata.rows) - fasta_genes):
        warn(f"gene {gene!r} in metadata TSV not found in FASTA: skipped")
    return merged


def _resolve_dbtype(flag: Dbtype | None, records: list[Record]) -> Dbtype:
    """Explicit --dbtype wins; else the abricate mol_type heuristic over the
    concatenated input sequences (db.mol_type, SPEC.md §2)."""
    if flag is not None:
        return flag
    return mol_type("".join(record.sequence for record in records))


def perform_build(
    name: str,
    fasta: Path,
    tsv: Path | None,
    dbtype: Dbtype | None,
    description: str,
    datadir: Path | None,
    force: bool,
    *,
    warn: Callable[[str], None],
    quiet: bool = True,
    kind: Kind | None = None,
    typing: Path | None = None,
    source: Literal["bundled"] | None = None,
) -> BuildReceipt:
    """Run the custom-build pipeline and return the receipt — the shared CLI
    + MCP path. Warnings go to the caller-supplied ``warn`` (CLI: stderr;
    MCP: dropped — stderr is reserved for the protocol).

    The input kind is detected by suffix (GBK/GFF -> cluster, else gene);
    an explicit ``kind`` must agree with the detection or the call fails as
    a usage error. Both branches take ``typing`` (a gapit.typing/1 or /2
    spec, validated — gene references against the FASTA records on the gene
    branch — then copied into the database); the FASTA-only ``tsv``/
    ``dbtype``/``description`` options are rejected on the cluster branch.
    """

    # Security/frozen rule: `Path(datadir) / name` REPLACES the base when name
    # is absolute (and `..` escapes it); plain names only, all else allowed.
    if not name or name in (".", "..") or "/" in name or "\\" in name:
        raise UsageError(
            f"database name must be a plain name without path separators: {name!r}",
            code="USAGE_ERROR",
            context={"name": name},
        )
    if not fasta.is_file():
        raise InputError(
            f"FASTA file not found or unreadable: {fasta}",
            code="INPUT_NOT_FOUND",
            context={"file": str(fasta)},
        )
    detected: Kind = "cluster" if is_cluster_input(fasta) else "gene"
    if kind is not None and kind != detected:
        raise UsageError(
            f"--kind {kind} contradicts the detected {detected} input format: {fasta.name}",
            code="USAGE_ERROR",
            context={"kind": kind, "detected": detected},
        )
    if detected == "cluster" and (tsv is not None or dbtype is not None or description):
        raise UsageError(
            "--tsv/--dbtype/--description apply to FASTA (gene) builds only",
            code="USAGE_ERROR",
        )
    db_dir = config.ensure_datadir(datadir) / name
    if (db_dir / "gapit-manifest.json").is_file() and not force:
        raise DatabaseError(
            f"won't overwrite existing database {name} (use --force)",
            code="DB_ALREADY_EXISTS",
            context={"db": name},
        )
    if detected == "cluster":
        manifest = perform_cluster_build(name, fasta, typing, db_dir, quiet=quiet)
        return BuildReceipt(
            db=name,
            records=manifest.n_records,
            dbtype=manifest.dbtype,
            destination=str(db_dir),
        )
    records = [_to_record(fasta_record, name, description) for fasta_record in iter_fasta(fasta)]
    if tsv is not None:
        records = _merge(records, _read_metadata(tsv, warn), warn)
    typing_schema = ""
    if typing is not None:
        document = read_typing_document(typing)
        validate_gene_typing(document, frozenset(record.gene for record in records))
        typing_schema = typing_schema_of(typing)
    db_dir.mkdir(parents=True, exist_ok=True)
    write_records(records, db_dir / "records.jsonl")
    if typing is not None:
        shutil.copyfile(typing, db_dir / "typing.json")
    manifest = build_database(
        db_dir,
        name=name,
        dbtype=_resolve_dbtype(dbtype, records),
        source_urls=("local",),
        fetched_at=datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        typing_schema=typing_schema,
        source=source,
        quiet=quiet,
    )
    return BuildReceipt(
        db=name,
        records=manifest.n_records,
        dbtype=manifest.dbtype,
        destination=str(db_dir),
    )
