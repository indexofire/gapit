"""The custom-database build path: `gapit db build` (from cmd_db.py).

Turns a user-supplied FASTA into a fully built gapit-native database via the
existing pipeline (records.jsonl -> sequences + BLAST index + .mmi + manifest,
written last). This module is orchestration plus a metadata merge only — the
building blocks live in records/dbbuild/fasta/dbcodec/db (SPEC.md §11 covers
the build pipeline itself).

Header kind is detected PER RECORD (mixed files allowed): a ``gapit|``
prefix decodes through the strict tagged codec, anything else through the
abricate ``~~~`` rules (a plain id carries no ``~~~`` and decodes to itself).
``db`` is always the TARGET name and ``source_id`` keeps the original id
token; record order is input order. Sequences are stored verbatim — no
provider-style normalization, this is the user's curated truth.
"""

import csv
from collections.abc import Callable
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Annotated, Literal

import typer
from pydantic import BaseModel

from gapit import config
from gapit.db import mol_type
from gapit.dbbuild import build_database
from gapit.dbcodec import decode_seqid
from gapit.errors import DatabaseError, GapitError, InputError, render_error
from gapit.fasta import FastaRecord, iter_fasta
from gapit.records import Record, write_records

Dbtype = Literal["nucl", "prot"]


class BuildReceipt(BaseModel, frozen=True):
    """One-line JSON success receipt for `db build` (mirrors db fetch)."""

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


def _dispatch(action: Callable[[], None]) -> None:
    """Run a command body; failures render the gapit.error/1 envelope on
    stderr and exit with the documented code (local mirror of cmd_db._dispatch
    — reportPrivateUsage blocks importing it, the Wave A3 ``_run`` precedent)."""
    try:
        action()
    except Exception as exc:
        typer.echo(render_error(exc), err=True)
        raise typer.Exit(code=exc.exit_code if isinstance(exc, GapitError) else 1) from exc


def _datadir_root(datadir: Path | None) -> Path:
    """Resolve the build datadir, creating it when absent (local mirror of
    cmd_db._fetch_root: `db build` is a write path like `db fetch`, so a
    fresh machine bootstraps instead of hitting DATADIR_NOT_FOUND)."""
    try:
        root = config.resolve_datadir(datadir)
    except DatabaseError as exc:
        if exc.code != "DATADIR_NOT_FOUND":
            raise
        root = Path(exc.context["datadir"])
    try:
        root.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise DatabaseError(
            f"cannot create datadir: {root}",
            code="DATADIR_CREATE_FAILED",
            context={"datadir": str(root)},
        ) from exc
    return root


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


def db_build_command(
    name: Annotated[
        str,
        typer.Argument(help="Target database name (created under the datadir)."),
    ],
    fasta: Annotated[
        Path,
        typer.Argument(
            help=(
                "Input FASTA: plain, abricate ~~~, or gapit| headers, detected per"
                " record (.gz/.bz2 accepted)."
            ),
        ),
    ],
    tsv: Annotated[
        Path | None,
        typer.Option(
            "--tsv", help="Metadata TSV: header row with gene/accession/function columns."
        ),
    ] = None,
    dbtype: Annotated[
        Dbtype | None,
        typer.Option("--dbtype", help="Force nucl or prot (default: auto-detect)."),
    ] = None,
    datadir: Annotated[
        Path | None,
        typer.Option(
            "--datadir",
            help="Database directory (default: $GAPIT_DATADIR, then ~/.local/share/gapit/db).",
        ),
    ] = None,
    description: Annotated[
        str,
        typer.Option(
            "--description",
            help="Default product for records whose FASTA header has no description text.",
        ),
    ] = "",
    force: Annotated[
        bool,
        typer.Option("--force", help="Overwrite the database if it already exists."),
    ] = False,
    quiet: Annotated[bool, typer.Option("--quiet", help="Silence stderr diagnostics.")] = False,
) -> None:
    """Build a custom gapit-native database from a FASTA (+ optional TSV)."""

    def warn(message: str) -> None:
        if not quiet:
            typer.echo(f"WARNING: {message}", err=True)

    def run() -> None:
        if not fasta.is_file():
            raise InputError(
                f"FASTA file not found or unreadable: {fasta}",
                code="INPUT_NOT_FOUND",
                context={"file": str(fasta)},
            )
        db_dir = _datadir_root(datadir) / name
        if (db_dir / "gapit-manifest.json").is_file() and not force:
            raise DatabaseError(
                f"won't overwrite existing database {name} (use --force)",
                code="DB_ALREADY_EXISTS",
                context={"db": name},
            )
        records = [
            _to_record(fasta_record, name, description) for fasta_record in iter_fasta(fasta)
        ]
        if tsv is not None:
            records = _merge(records, _read_metadata(tsv, warn), warn)
        db_dir.mkdir(parents=True, exist_ok=True)
        write_records(records, db_dir / "records.jsonl")
        manifest = build_database(
            db_dir,
            name=name,
            dbtype=_resolve_dbtype(dbtype, records),
            source_urls=("local",),
            fetched_at=datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
            quiet=quiet,
        )
        typer.echo(
            BuildReceipt(
                db=name,
                records=manifest.n_records,
                dbtype=manifest.dbtype,
                destination=str(db_dir),
            ).model_dump_json()
        )

    _dispatch(run)
