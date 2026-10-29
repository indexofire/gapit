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
- ``db info NAME``: one database's provenance, literature references, intended
  usage, and interpretation caveats (use-case: :mod:`gapit.db_info_ops`,
  command: :mod:`gapit.cmd_db_info`).
- ``db list``: show every known database with its NAME, DESCRIPTION, and
  installed state (PROVIDER last; rich table on a TTY, byte-stable TSV
  otherwise; use-case: :mod:`gapit.db_ops`).
- ``db outdated`` / ``db search``: read-only queries over the installed
  databases (use-case: :mod:`gapit.db_query_ops`; commands:
  :mod:`gapit.cmd_db_outdated`, :mod:`gapit.cmd_db_search`).
- ``db rename OLD NEW``: re-home an installed database under a new name
  (use-case: :mod:`gapit.db_rename_ops`, command:
  :mod:`gapit.cmd_db_rename`) — the migration path for names that predate
  the naming convention.
- ``db install``: checksum-verified LOCAL-FILE installation
  (:mod:`gapit.cmd_db_install`, split off in Wave G at the 250 LOC ceiling).
"""

import sys
from pathlib import Path

import click
from rich.console import Console

from gapit import config
from gapit.cmd_db_build import db_build_command
from gapit.cmd_db_info import db_info_command
from gapit.cmd_db_install import db_install_command
from gapit.cmd_db_outdated import db_outdated_command
from gapit.cmd_db_rename import db_rename_command
from gapit.cmd_db_search import db_search_command
from gapit.db_ops import (
    datadir_json,
    db_list_entries,
    db_list_json,
    db_list_status,
    db_list_table,
    db_list_typing,
    perform_fetch,
)
from gapit.dispatch import datadir_option, dispatch, quiet_option
from gapit.errors import usage_fail
from gapit.formats.csvutil import tsv_cell


@click.argument(
    "name",
    required=False,
    help=("Database name (see: gapit db list), or 'all' for the default set (card, vfdb)."),
)
@datadir_option()
@click.option("--force", "-f", is_flag=True, help="Overwrite the database if it already exists.")
@quiet_option()
@click.option("--debug", "-v", is_flag=True, help="Echo external command lines to stderr.")
def db_fetch_command(
    name: str | None,
    datadir: Path | None,
    force: bool,
    quiet: bool,
    debug: bool,
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
        for _receipt in perform_fetch(
            name,
            datadir,
            force=force,
            quiet=quiet,
            debug=debug,
            emit=lambda receipt: click.echo(receipt.model_dump_json()),
        ):
            pass

    dispatch(run)


@datadir_option()
@click.option(
    "--json", "-J", "as_json", is_flag=True, help="Print the gapit.datadir/1 JSON document."
)
def db_path_command(datadir: Path | None, as_json: bool) -> None:
    """Print the resolved datadir (the directory `db list`/`db fetch` use)."""

    def run() -> None:
        if as_json:
            click.echo(datadir_json(datadir))
            return
        click.echo(str(config.datadir_path(datadir)))

    dispatch(run)


@datadir_option()
@click.option(
    "--json", "-J", "as_json", is_flag=True, help="Print machine-readable JSON instead of a table."
)
def db_list_command(datadir: Path | None, as_json: bool) -> None:
    """List known databases and their state.

    Columns: NAME, DESCRIPTION, STATUS, DBTYPE, TYPING, PROVIDER.
    """

    def run() -> None:
        # Catalog bootstrap: a fresh pip/conda install has no datadir yet, but
        # listing must still work and leave the default location ready for the
        # first fetch/screen. Installed-state checks below only read manifests.
        root = config.ensure_datadir(datadir)
        if as_json:
            click.echo(db_list_json(root))
            return
        entries = db_list_entries(root)
        if sys.stdout.isatty():
            # TTY-only rich table; pipes/redirects must keep the byte-stable TSV.
            Console().print(db_list_table(entries))
            return
        click.echo("NAME\tDESCRIPTION\tSTATUS\tDBTYPE\tTYPING\tPROVIDER")
        for entry in entries:
            click.echo(
                f"{tsv_cell(entry.name)}\t{tsv_cell(entry.description)}"
                f"\t{db_list_status(entry)}\t{entry.dbtype}"
                f"\t{db_list_typing(entry)}\t{tsv_cell(entry.vendor)}"
            )

    dispatch(run)


def register_db_command(app: click.Group) -> None:
    """Attach the `db` command group (install, fetch, list) to the CLI app.

    Input-requiring subcommands print help when invoked bare
    (``no_args_is_help``); ``list``/``outdated`` are valid with no
    arguments and keep their behavior.
    """
    db_app = click.Group(
        "db",
        help=(
            "Database acquisition and maintenance: fetch, build, verified local install,"
            " listing, per-database info, rename, search, staleness."
        ),
        no_args_is_help=True,
    )
    db_app.command("install", no_args_is_help=True)(db_install_command)
    db_app.command("build", no_args_is_help=True)(db_build_command)
    db_app.command("fetch", no_args_is_help=True)(db_fetch_command)
    db_app.command("path")(db_path_command)
    db_app.command("list")(db_list_command)
    db_app.command("info", no_args_is_help=True)(db_info_command)
    db_app.command("outdated")(db_outdated_command)
    db_app.command("search", no_args_is_help=True)(db_search_command)
    db_app.command("rename", no_args_is_help=True)(db_rename_command)
    app.add_command(db_app, "db")
