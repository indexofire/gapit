"""The `gapit db outdated` command: installed-database staleness report.

Reports every installed database's age against a staleness threshold and
the bundled snapshot date. Staleness is a REPORT, never an error state:
exit 0 even when everything is stale. Read-only: no builds, no network,
no datadir writes.
"""

from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path
from typing import Annotated, Literal

import typer
from pydantic import BaseModel, ConfigDict, Field

from gapit import config
from gapit.dispatch import Datadir, dispatch
from gapit.errors import DatabaseError, InputError
from gapit.providers import REGISTRY
from gapit.providers.common import bundled_snapshot_manifest
from gapit.records import installed_db_dirs, read_manifest

DEFAULT_STALE_DAYS = 90


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


def db_outdated_command(
    datadir: Datadir = None,
    days: Annotated[
        int,
        typer.Option("--days", min=0, help="Staleness threshold in days (stale past this age)."),
    ] = DEFAULT_STALE_DAYS,
    as_json: Annotated[
        bool,
        typer.Option("--json", help="Print machine-readable JSON instead of a table."),
    ] = False,
    quiet: Annotated[bool, typer.Option("--quiet", help="Silence stderr diagnostics.")] = False,
) -> None:
    """Report installed database ages and available updates (exit 0 however
    stale things are — a report, not an error).

    A database is `stale` past --days (default 90) and `snapshot-update`
    when its provider's bundled snapshot is newer than the installed copy.
    """

    def run() -> None:
        entries = perform_outdated(datadir, days=days)
        if as_json:
            document = DbOutdatedDocument(databases=tuple(entries))
            typer.echo(document.model_dump_json(indent=2, by_alias=True))
            return
        for line in outdated_tsv_lines(entries):
            typer.echo(line)

    dispatch(run)
