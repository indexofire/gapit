"""The `gapit db search` command: typer shell over the use-case
(:mod:`gapit.db_query_ops`) — case-insensitive lookup across the
``records.jsonl`` truth stores of every installed database (SPEC.md §11).

Read-only: no builds, no network, no datadir writes.
"""

from typing import Annotated

import typer

from gapit.db_query_ops import (
    DEFAULT_LIMIT,
    SearchField,
    perform_search,
    search_json_line,
    search_tsv_row,
)
from gapit.dispatch import Datadir, dispatch
from gapit.proctools import note


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
            render=search_json_line if as_json else search_tsv_row,
            quiet=quiet,
        )
        if limit > 0 and total > limit:
            note(quiet, f"truncated to {limit} of {total} matching records (use --limit 0 for all)")
        for line in hits:
            typer.echo(line)

    dispatch(run)
