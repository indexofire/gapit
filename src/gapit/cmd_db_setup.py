"""The `gapit db setup` command: bootstrap the datadir and rebuild indices.

With no selection, creates a missing datadir, materializes bundled databases
that are not installed, then rebuilds the BLAST index for every database under
the datadir. ``--db/-d NAME`` restricts the work to named databases (a missing
provider is fetched, a missing bundle materialized); ``--all/-a`` explicitly
selects every installed/bundled database. ``--json`` emits the
``gapit.dbsetup/1`` receipt; ``--no-bundles`` skips bundle materialization.
"""

from pathlib import Path

import click

from gapit.db_setup_ops import perform_setup
from gapit.dispatch import datadir_option, dispatch, quiet_option
from gapit.errors import usage_fail


@click.command("setup")
@click.option(
    "--db",
    "-d",
    "names",
    multiple=True,
    metavar="NAME",
    help=(
        "Database to set up (repeatable). A missing provider is fetched; a"
        " missing bundled database is materialized."
    ),
)
@click.option(
    "--all",
    "-a",
    "all_databases",
    is_flag=True,
    help="Set up every installed/bundled database (never fetches providers).",
)
@datadir_option()
@quiet_option()
@click.option("--debug", "-v", is_flag=True, help="Echo external command lines to stderr.")
@click.option(
    "--no-bundles",
    "-B",
    is_flag=True,
    help="Skip bundled materialization; named providers may still be fetched.",
)
@click.option(
    "--json",
    "-J",
    "as_json",
    is_flag=True,
    help="Print the gapit.dbsetup/1 JSON document (progress is suppressed).",
)
def db_setup_command(
    names: tuple[str, ...],
    all_databases: bool,
    datadir: Path | None,
    quiet: bool,
    debug: bool,
    no_bundles: bool,
    as_json: bool,
) -> None:
    """Set up all databases, or only the selected --db NAME(s), under the datadir.

    With neither --all nor --db, behaves like --all. --db may fetch a missing
    provider and materialize a missing bundled database; --all never fetches.
    """

    def run() -> None:
        if all_databases and names:
            usage_fail("db setup takes --all or --db NAME, not both")
        # Keep the bare invocation as the historical all-databases path.
        selected_all = all_databases or not names
        # --json is an agent-facing receipt: keep stderr clean as well, so
        # the caller sees exactly one JSON document and no progress noise.
        document = perform_setup(
            datadir,
            debug=debug,
            quiet=quiet or as_json,
            bundles=not no_bundles,
            names=None if selected_all else names,
        )
        if as_json:
            click.echo(document.model_dump_json(by_alias=True))

    dispatch(run)
