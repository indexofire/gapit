"""The `gapit db` command group: build, fetch, list, install.

- ``db build NAME FASTA``: build a custom gapit-native database from a
  user-supplied FASTA (+ optional TSV metadata) — the whole acquisition
  pipeline without a provider (use-case: :mod:`gapit.db_build_ops`,
  command: :mod:`gapit.cmd_db_build`).
- ``db fetch NAME|all``: acquire database(s) under the datadir. The literal
  NAME ``all`` installs the default set (``DEFAULT_DBS``: card, vfdb) in
  order — every provider downloads from its upstream source at fetch time
  (nothing is bundled: several upstream licenses forbid redistribution).
  One JSON receipt line per database on stdout (use-case:
  :mod:`gapit.db_ops`). Bare ``db fetch`` prints help; a download is never
  implied by an omitted NAME.
- ``db list``: show every known database with its NAME, upstream PROVIDER,
  and installed state (rich table on a TTY, byte-stable TSV otherwise;
  use-case: :mod:`gapit.db_ops`).
- ``db outdated`` / ``db search``: read-only queries over the installed
  databases (use-case: :mod:`gapit.db_query_ops`; commands:
  :mod:`gapit.cmd_db_outdated`, :mod:`gapit.cmd_db_search`).
- ``db install``: checksum-verified LOCAL-FILE installation
  (:mod:`gapit.cmd_db_install`, split off in Wave G at the 250 LOC ceiling).
"""

import sys
from typing import Annotated

import typer
from rich.console import Console

from gapit import config
from gapit.cmd_db_build import db_build_command
from gapit.cmd_db_install import db_install_command
from gapit.cmd_db_outdated import db_outdated_command
from gapit.cmd_db_search import db_search_command
from gapit.db_ops import db_list_entries, db_list_json, db_list_status, db_list_table, perform_fetch
from gapit.dispatch import Datadir, dispatch
from gapit.errors import usage_fail


def db_fetch_command(
    name: Annotated[
        str | None,
        typer.Argument(
            help=("Database name (see: gapit db list), or 'all' for the default set (card, vfdb)."),
        ),
    ] = None,
    datadir: Datadir = None,
    force: Annotated[
        bool,
        typer.Option("--force", help="Overwrite the database if it already exists."),
    ] = False,
    quiet: Annotated[bool, typer.Option("--quiet", help="Silence stderr diagnostics.")] = False,
    debug: Annotated[
        bool,
        typer.Option("--debug", help="Echo external command lines to stderr."),
    ] = False,
) -> None:
    """Fetch and build database(s) into <datadir>/NAME.

    NAME 'all' installs every database in DEFAULT_DBS order; each success
    prints its own one-line JSON receipt to stdout (stdout purity: receipts
    are data).
    """

    def run() -> None:
        # Bare `db fetch` prints help (no_args_is_help); only a flags-present
        # invocation with no NAME reaches this guard.
        if name is None:
            usage_fail("db fetch requires a database NAME or 'all' (see: gapit db list)")
        for receipt in perform_fetch(name, datadir, force=force, quiet=quiet, debug=debug):
            typer.echo(receipt.model_dump_json())

    dispatch(run)


def db_list_command(
    datadir: Datadir = None,
    as_json: Annotated[
        bool,
        typer.Option("--json", help="Print machine-readable JSON instead of a table."),
    ] = False,
) -> None:
    """List known databases (NAME, upstream PROVIDER) and their installed state."""

    def run() -> None:
        root = config.resolve_datadir(datadir)
        if as_json:
            typer.echo(db_list_json(root))
            return
        entries = db_list_entries(root)
        if sys.stdout.isatty():
            # TTY-only rich table; pipes/redirects must keep the byte-stable TSV.
            Console().print(db_list_table(entries))
            return
        typer.echo("NAME\tPROVIDER\tSTATUS\tDBTYPE\tDESCRIPTION")
        for entry in entries:
            typer.echo(
                f"{entry.name}\t{entry.vendor}\t{db_list_status(entry)}"
                f"\t{entry.dbtype}\t{entry.description}"
            )

    dispatch(run)


def register_db_command(app: typer.Typer) -> None:
    """Attach the `db` command group (install, fetch, list) to the CLI app.

    Input-requiring subcommands print help when invoked bare
    (``no_args_is_help``); ``list``/``outdated`` are valid with no
    arguments and keep their behavior.
    """
    db_app = typer.Typer(
        help="Database acquisition and maintenance (database fetch, verified local install).",
        no_args_is_help=True,
    )
    db_app.command("install", no_args_is_help=True)(db_install_command)
    db_app.command("build", no_args_is_help=True)(db_build_command)
    db_app.command("fetch", no_args_is_help=True)(db_fetch_command)
    db_app.command("list")(db_list_command)
    db_app.command("outdated")(db_outdated_command)
    db_app.command("search", no_args_is_help=True)(db_search_command)
    app.add_typer(db_app, name="db")
