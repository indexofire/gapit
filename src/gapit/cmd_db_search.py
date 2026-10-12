"""The `gapit db search` command: click shell over the use-case
(:mod:`gapit.db_query_ops`) — case-insensitive lookup across the
``records.jsonl`` truth stores of every installed database (SPEC.md §11).

Read-only: no builds, no network, no datadir writes.
"""

from pathlib import Path

import click

from gapit.db_query_ops import (
    DEFAULT_LIMIT,
    SearchField,
    perform_search,
    search_json_line,
    search_tsv_row,
)
from gapit.dispatch import EnumChoice, datadir_option, dispatch, quiet_option
from gapit.proctools import note


@click.argument("term", help="Search term (case-insensitive).")
@datadir_option()
@click.option(
    "--db",
    "-d",
    "-db",
    help="Restrict the scan to one installed database.",
)
@click.option(
    "--field",
    "-f",
    type=EnumChoice(SearchField),
    default=SearchField.any,
    help="Field to match: gene, accession, function, product, any.",
)
@click.option("--exact", "-e", is_flag=True, help="Full-field equality instead of substring.")
@click.option(
    "--limit",
    "-l",
    type=click.IntRange(min=0),
    default=DEFAULT_LIMIT,
    help="Max hits to print (0 = unlimited).",
)
@click.option(
    "--json",
    "-J",
    "as_json",
    is_flag=True,
    help="Print JSONL lines instead of TSV rows.",
)
@quiet_option()
def db_search_command(
    term: str,
    datadir: Path | None,
    db: str | None,
    field: SearchField,
    exact: bool,
    limit: int,
    as_json: bool,
    quiet: bool,
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
            render=search_json_line if as_json else search_tsv_row,
            quiet=quiet,
        )
        if limit > 0 and total > limit:
            note(quiet, f"truncated to {limit} of {total} matching records (use --limit 0 for all)")
        for line in hits:
            click.echo(line)

    dispatch(run)
