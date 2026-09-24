"""The `gapit db outdated` command: typer shell over the use-case
(:mod:`gapit.db_query_ops`) — installed-database staleness report.

Reports every installed database's age against a staleness threshold and
the bundled snapshot date. Staleness is a REPORT, never an error state:
exit 0 even when everything is stale. Read-only: no builds, no network,
no datadir writes.
"""

from typing import Annotated

import typer

from gapit.db_query_ops import (
    DEFAULT_STALE_DAYS,
    DbOutdatedDocument,
    outdated_tsv_lines,
    perform_outdated,
)
from gapit.dispatch import Datadir, dispatch


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
