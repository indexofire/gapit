"""Read-only queries over installed databases, shared by the CLI and MCP.

- ``perform_search``: case-insensitive lookup across the ``records.jsonl``
  truth stores of every installed database (SPEC.md §11).
- ``perform_outdated``: installed-database staleness report against a
  threshold (``days``) and the bundled snapshot dates. Staleness is a
  REPORT, never an error state.

Both are read-only: no builds, no network, no datadir writes. The typer
commands (:mod:`gapit.cmd_db_search`, :mod:`gapit.cmd_db_outdated`) and the
MCP tools (:mod:`gapit.mcp_tools`) are thin callers.
"""

import enum
from collections.abc import Callable, Iterable
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from gapit import config
from gapit.errors import DatabaseError, InputError, UsageError
from gapit.proctools import note
from gapit.providers import REGISTRY
from gapit.providers.common import bundled_snapshot_manifest
from gapit.records import Record, installed_db_dirs, read_manifest, read_records

DEFAULT_LIMIT = 100
DEFAULT_STALE_DAYS = 90


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


def search_tsv_row(record: Record) -> str:
    """The hit row: DB\tGENE\tACCESSION\tFUNCTION\tPRODUCT\tLENGTH."""
    return (
        f"{record.db}\t{record.gene}\t{record.accession}\t"
        f"{';'.join(record.function)}\t{record.product}\t{len(record.sequence)}"
    )


def search_json_line(record: Record) -> str:
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
    render: Callable[[Record], str] = search_tsv_row,
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


class DbOutdatedEntry(BaseModel, frozen=True):
    """One installed-database row in the gapit.dboutdated/1 report."""

    db: str
    fetched_at: str
    age_days: float
    status: Literal["ok", "stale", "snapshot-update", "stale+snapshot-update"]
    upstream_version: str


class DbOutdatedDocument(BaseModel, frozen=True):
    """gapit.dboutdated/1 — `gapit db outdated --json` output.

    A CLI listing, deliberately local to this module: NOT registered in
    `gapit schema` (mirrors the gapit.dblist/1 precedent).
    """

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.dboutdated/1"] = Field(default="gapit.dboutdated/1", alias="schema")
    databases: tuple[DbOutdatedEntry, ...]


def _utc_timestamp(value: str, where: str) -> datetime:
    """Parse a manifest fetched_at into aware UTC; anything else is a
    malformed manifest (trusted metadata: writers emit ISO-8601 with Z)."""
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError as exc:
        raise InputError(
            f"malformed manifest {where}: invalid fetched_at: {value}",
            code="MANIFEST_MALFORMED",
            context={"file": where},
        ) from exc
    if parsed.tzinfo is None:
        raise InputError(
            f"malformed manifest {where}: fetched_at has no UTC offset: {value}",
            code="MANIFEST_MALFORMED",
            context={"file": where},
        )
    return parsed.astimezone(UTC)


def _snapshot_newer(name: str, fetched: datetime) -> bool:
    """True when the provider's bundled snapshot manifest is newer than the
    installed ``fetched`` (unknown provider / no bundle -> False)."""
    provider = REGISTRY.get(name)
    if provider is None:
        return False
    archived = bundled_snapshot_manifest(provider)
    return archived is not None and fetched < _utc_timestamp(
        archived.fetched_at, f"bundled snapshot of {name}"
    )


def perform_outdated(
    datadir: Path | None, *, days: int = DEFAULT_STALE_DAYS
) -> list[DbOutdatedEntry]:
    """Compute the staleness report entries — the shared CLI + MCP path.

    A database is `stale` past ``days`` and `snapshot-update` when its
    provider's bundled snapshot is newer than the installed copy.
    """
    root = config.resolve_datadir(datadir)
    db_dirs = installed_db_dirs(root)
    if not db_dirs:
        raise DatabaseError(
            f"no installed databases in datadir: {root}",
            code="DATADIR_EMPTY",
            context={"datadir": str(root)},
        )
    now = datetime.now(UTC)
    entries: list[DbOutdatedEntry] = []
    for db_dir in db_dirs:
        manifest_path = db_dir / "gapit-manifest.json"
        manifest = read_manifest(manifest_path)
        fetched = _utc_timestamp(manifest.fetched_at, str(manifest_path))
        age = (now - fetched).total_seconds() / 86400
        match (age > days, _snapshot_newer(db_dir.name, fetched)):
            case (True, True):
                status = "stale+snapshot-update"
            case (True, False):
                status = "stale"
            case (False, True):
                status = "snapshot-update"
            case (False, False):
                status = "ok"
        entries.append(
            DbOutdatedEntry(
                db=db_dir.name,
                fetched_at=manifest.fetched_at,
                age_days=round(age, 2),
                status=status,
                upstream_version=manifest.upstream_version,
            )
        )
    return entries


def outdated_tsv_lines(entries: Iterable[DbOutdatedEntry]) -> list[str]:
    """The NAME/FETCHED_AT/AGE_DAYS/STATUS table lines (CLI stdout = MCP text)."""
    lines = ["NAME\tFETCHED_AT\tAGE_DAYS\tSTATUS"]
    lines.extend(
        f"{entry.db}\t{entry.fetched_at}\t{entry.age_days:.2f}\t{entry.status}" for entry in entries
    )
    return lines
