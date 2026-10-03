"""The custom-database build use-case, shared by the CLI and MCP.

Turns a user-supplied input into a fully built gapit-native database: FASTA
inputs run the gene pipeline (records.jsonl -> sequences + BLAST index +
manifest, written last), GBK/GFF inputs run the cluster pipeline (locus
FASTA + features.json + manifest ``kind: cluster``, gapit.dbbuild). This
module is orchestration plus kind dispatch only — the building blocks live in
records/dbbuild/db_build_meta/fasta/dbcodec/db (SPEC.md §11). The
typer command (:mod:`gapit.cmd_db_build`) and the MCP tool
(:mod:`gapit.mcp_tools`) are thin callers.

Header kind is detected PER RECORD (mixed files allowed): a ``gapit|``
prefix decodes through the strict tagged codec, anything else through the
abricate ``~~~`` rules (a plain id carries no ``~~~`` and decodes to itself).
``db`` is always the TARGET name and ``source_id`` keeps the original id
token; record order is input order. Sequences are stored verbatim — no
provider-style normalization, this is the user's curated truth.
"""

import shutil
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

from pydantic import BaseModel

from gapit import config
from gapit.clusterbuild import is_cluster_input, perform_cluster_build
from gapit.db import mol_type
from gapit.db_build_meta import merge_metadata, read_metadata, split_function
from gapit.dbbuild import build_database, drop_exact_duplicates
from gapit.dbcodec import decode_seqid
from gapit.errors import DatabaseError, InputError, UsageError
from gapit.fasta import FastaRecord, iter_fasta
from gapit.gene_floors import FLOORS_FILENAME, read_floors, validate_floors
from gapit.proctools import note
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
    # Exact (gene, sequence) duplicates dropped at build time; None (omitted
    # at serialization) when none were, so no-dup receipts stay
    # byte-identical to the pre-field shape.
    duplicates_dropped: int | None = None


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
        function=split_function(header.function),
        product=fasta.description or default_product or header.gene,
        source_id=fasta.id,
    )


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
    floors: Path | None = None,
    source: Literal["bundled"] | None = None,
) -> BuildReceipt:
    """Run the custom-build pipeline and return the receipt — the shared CLI
    + MCP path. Warnings go to the caller-supplied ``warn`` (CLI: stderr;
    MCP: dropped — stderr is reserved for the protocol).

    The input kind is detected by suffix (GBK/GFF -> cluster, else gene);
    an explicit ``kind`` must agree with the detection or the call fails as
    a usage error. Both branches take ``typing`` (a gapit.typing/1 or /2
    spec, validated — gene references against the FASTA records on the gene
    branch — then copied into the database); the gene branch also takes
    ``floors`` (a gapit.floors/1 per-gene identity-floor document, validated
    the same way and copied in as ``floors.json`` — the file's presence is
    the flag, no manifest field); the FASTA-only ``tsv``/``dbtype``/
    ``description``/``floors`` options are rejected on the cluster branch.

    The gene branch drops exact ``(gene, sequence)`` duplicates after
    parsing (first kept; see :func:`gapit.dbbuild.drop_exact_duplicates`)
    — the count lands in the receipt and a note on stderr when not quiet.
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
    if detected == "cluster" and floors is not None:
        raise UsageError(
            "--floors apply to FASTA (gene) builds only",
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
        records = merge_metadata(records, read_metadata(tsv, warn), warn)
    records, duplicates_dropped = drop_exact_duplicates(records)
    if duplicates_dropped:
        note(
            quiet,
            f"dropped {duplicates_dropped} exact duplicate record(s) "
            "(same gene+sequence; first kept)",
        )
    typing_schema = ""
    if typing is not None:
        document = read_typing_document(typing)
        validate_gene_typing(document, frozenset(record.gene for record in records))
        typing_schema = typing_schema_of(typing)
    if floors is not None:
        validate_floors(
            read_floors(floors), frozenset(record.gene for record in records), source=floors
        )
    db_dir.mkdir(parents=True, exist_ok=True)
    write_records(records, db_dir / "records.jsonl")
    if typing is not None:
        shutil.copyfile(typing, db_dir / "typing.json")
    if floors is not None:
        shutil.copyfile(floors, db_dir / FLOORS_FILENAME)
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
        duplicates_dropped=duplicates_dropped or None,
    )
