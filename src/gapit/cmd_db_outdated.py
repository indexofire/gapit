"""The `gapit db outdated` command: click shell over the use-case
(:mod:`gapit.db_query_ops`) — installed-database staleness report.

Reports every installed database's age against a staleness threshold.
Staleness is a REPORT, never an error state: exit 0 even when everything is
stale. Read-only: no builds, no network, no datadir writes.
"""

from pathlib import Path

import click

from gapit.db_query_ops import (
    DEFAULT_STALE_DAYS,
    DbOutdatedDocument,
    outdated_tsv_lines,
    perform_outdated,
)
from gapit.dispatch import datadir_option, dispatch, quiet_option


@datadir_option()
@click.option(
    "--days",
    "-d",
    type=click.IntRange(min=0),
    default=DEFAULT_STALE_DAYS,
    help="Staleness threshold in days (stale past this age).",
)
@click.option(
    "--json", "-J", "as_json", is_flag=True, help="Print machine-readable JSON instead of a table."
)
@quiet_option()
def db_outdated_command(datadir: Path | None, days: int, as_json: bool, quiet: bool) -> None:
    """Report installed database ages and available updates (exit 0 however
    stale things are — a report, not an error).

    A database is `stale` past --days (default 90), `ok` otherwise.
    """

    def run() -> None:
        entries = perform_outdated(datadir, days=days)
        if as_json:
            document = DbOutdatedDocument(databases=tuple(entries))
            click.echo(document.model_dump_json(indent=2, by_alias=True))
            return
        for line in outdated_tsv_lines(entries):
            click.echo(line)

    dispatch(run)
