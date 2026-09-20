"""The `gapit db search` command: case-insensitive lookup across the
``records.jsonl`` truth stores of every installed database (SPEC.md §11).

Read-only: no builds, no network, no datadir writes.
"""

import enum
from collections.abc import Callable
from pathlib import Path
from typing import Annotated

import typer
from pydantic import BaseModel

from gapit import config
from gapit.dispatch import Datadir, dispatch
from gapit.errors import DatabaseError, UsageError
from gapit.proctools import note
from gapit.records import Record, installed_db_dirs, read_records

DEFAULT_LIMIT = 100


class SearchField(enum.Enum):
    """Record fields `db search` can match against."""

    gene = "gene"
    accession = "accession"
    function = "function"
    product = "product"
    any = "any"


class SearchEntry(BaseModel, frozen=True):
    """One search hit: the JSONL line (and the TSV row's source), a CLI
    listing shape — not a registered schema."""

    db: str
    gene: str
    accession: str
    function: tuple[str, ...]
    product: str
    length: int


def _candidates(record: Record, field: SearchField) -> tuple[str, ...]:
    """The string values of the chosen field (function contributes one per class)."""
    match field:
        case SearchField.gene:
            return (record.gene,)
        case SearchField.accession:
            return (record.accession,)
        case SearchField.function:
            return record.function
        case SearchField.product:
            return (record.product,)
        case SearchField.any:
            return (record.gene, record.accession, *record.function, record.product)


def _matches(record: Record, needle: str, field: SearchField, exact: bool) -> bool:
    """Case-insensitive substring by default; --exact is full-field equality."""
    values = _candidates(record, field)
    if exact:
        return any(value.lower() == needle for value in values)
    return any(needle in value.lower() for value in values)


def _tsv_row(record: Record) -> str:
    """The hit row: DB\tGENE\tACCESSION\tFUNCTION\tPRODUCT\tLENGTH."""
    return (
        f"{record.db}\t{record.gene}\t{record.accession}\t"
        f"{';'.join(record.function)}\t{record.product}\t{len(record.sequence)}"
    )


def _json_line(record: Record) -> str:
    """The hit as one JSONL line (same fields as the TSV row, snake_case)."""
    return SearchEntry(
        db=record.db,
        gene=record.gene,
        accession=record.accession,
        function=record.function,
        product=record.product,
        length=len(record.sequence),
    ).model_dump_json()


def perform_search(
    term: str,
    datadir: Path | None,
    *,
    db: str | None = None,
    field: SearchField = SearchField.any,
    exact: bool = False,
    limit: int = DEFAULT_LIMIT,
    render: Callable[[Record], str] = _tsv_row,
    quiet: bool = True,
) -> tuple[list[str], int]:
    """Scan the installed records.jsonl stores — the shared CLI + MCP path.

    Returns ``(hit lines up to limit, total matching count)``; the caller
    owns the truncation note (a stderr diagnostic) and the printing.
    """
    root = config.resolve_datadir(datadir)
    db_dirs = installed_db_dirs(root)
    if db is not None and db not in {db_dir.name for db_dir in db_dirs}:
        installed = " ".join(sorted(db_dir.name for db_dir in db_dirs))
        raise UsageError(
            f"unknown database: {db} (installed: {installed or 'none'})",
            code="USAGE_ERROR",
            context={"db": db, "installed": installed},
        )
    needle = term.lower()
    hits: list[str] = []
    total = 0
    for db_dir in db_dirs:
        if db is not None and db_dir.name != db:
            continue
        records_path = db_dir / "records.jsonl"
        if not records_path.is_file():
            if db is not None:
                raise DatabaseError(
                    f"database {db_dir.name} is incomplete: records.jsonl missing",
                    code="DB_INCOMPLETE",
                    context={"db": db_dir.name},
                )
            note(quiet, f"skipping {db_dir.name}: no records.jsonl (incomplete database)")
            continue
        for record in read_records(records_path):
            if _matches(record, needle, field, exact):
                total += 1
                if limit == 0 or len(hits) < limit:
                    hits.append(render(record))
    return hits, total


def db_search_command(
    term: Annotated[str, typer.Argument(help="Search term (case-insensitive).")],
    datadir: Datadir = None,
    db: Annotated[
        str | None,
        typer.Option("--db", help="Restrict the scan to one installed database."),
    ] = None,
    field: Annotated[
        SearchField,
        typer.Option("--field", help="Field to match: gene, accession, function, product, any."),
    ] = SearchField.any,
    exact: Annotated[
        bool,
        typer.Option("--exact", help="Full-field equality instead of substring."),
    ] = False,
    limit: Annotated[
        int,
        typer.Option("--limit", min=0, help="Max hits to print (0 = unlimited)."),
    ] = DEFAULT_LIMIT,
    as_json: Annotated[
        bool,
        typer.Option("--json", help="Print JSONL lines instead of TSV rows."),
    ] = False,
    quiet: Annotated[bool, typer.Option("--quiet", help="Silence stderr diagnostics.")] = False,
) -> None:
    """Search records.jsonl across installed databases (zero hits exit 0)."""

    def run() -> None:
        hits, total = perform_search(
            term,
            datadir,
            db=db,
            field=field,
            exact=exact,
            limit=limit,
            render=_json_line if as_json else _tsv_row,
            quiet=quiet,
        )
        if limit > 0 and total > limit:
            note(quiet, f"truncated to {limit} of {total} matching records (use --limit 0 for all)")
        for line in hits:
            typer.echo(line)

    dispatch(run)
