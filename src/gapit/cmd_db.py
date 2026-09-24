"""The `gapit db` command group: build, fetch, list, install.

- ``db build NAME FASTA``: build a custom gapit-native database from a
  user-supplied FASTA (+ optional TSV metadata) — the whole acquisition
  pipeline without a provider (use-case: :mod:`gapit.db_build_ops`,
  command: :mod:`gapit.cmd_db_build`).
- ``db fetch [NAME]``: acquire provider database(s) under the datadir. With
  no NAME the default set (``DEFAULT_DBS``: card, vfdb) installs in order —
  each from its BUNDLED SNAPSHOT (Wave G) when one resolves, else over the
  network; ``--from-source`` forces the upstream download even when a
  snapshot exists. One JSON receipt line per database on stdout
  (use-case: :mod:`gapit.db_ops`).
- ``db list``: show every known provider with its installed state
  (use-case: :mod:`gapit.db_ops`).
- ``db outdated`` / ``db search``: read-only queries over the installed
  databases (use-case: :mod:`gapit.db_query_ops`; commands:
  :mod:`gapit.cmd_db_outdated`, :mod:`gapit.cmd_db_search`).
- ``db install``: checksum-verified LOCAL-FILE installation
  (:mod:`gapit.cmd_db_install`, split off in Wave G at the 250 LOC ceiling).
"""

from typing import Annotated

import typer

from gapit import config
from gapit.cmd_db_build import db_build_command
from gapit.cmd_db_install import db_install_command
from gapit.cmd_db_outdated import db_outdated_command
from gapit.cmd_db_search import db_search_command
from gapit.db_ops import db_list_entries, db_list_json, perform_fetch
from gapit.dispatch import Datadir, dispatch


def db_fetch_command(
    name: Annotated[
        str | None,
        typer.Argument(
            help=(
                "Provider name (see: gapit db list). Omitted: install the default set "
                "(card, vfdb) from their bundled snapshots."
            ),
        ),
    ] = None,
    datadir: Datadir = None,
    force: Annotated[
        bool,
        typer.Option("--force", help="Overwrite the database if it already exists."),
    ] = False,
    from_source: Annotated[
        bool,
        typer.Option(
            "--from-source", help="Ignore bundled snapshots, download from upstream sources."
        ),
    ] = False,
    quiet: Annotated[bool, typer.Option("--quiet", help="Silence stderr diagnostics.")] = False,
    debug: Annotated[
        bool,
        typer.Option("--debug", help="Echo external command lines to stderr."),
    ] = False,
) -> None:
    """Fetch and build provider database(s) into <datadir>/NAME.

    NAME omitted installs every database in DEFAULT_DBS order; each success
    prints its own one-line JSON receipt to stdout (stdout purity: receipts
    are data).
    """

    def run() -> None:
        for receipt in perform_fetch(
            name, datadir, force=force, from_source=from_source, quiet=quiet, debug=debug
        ):
            typer.echo(receipt.model_dump_json())

    dispatch(run)


def db_list_command(
    datadir: Datadir = None,
    as_json: Annotated[
        bool,
        typer.Option("--json", help="Print machine-readable JSON instead of a table."),
    ] = False,
) -> None:
    """List database providers and their installed state under the datadir."""

    def run() -> None:
        root = config.resolve_datadir(datadir)
        if as_json:
            typer.echo(db_list_json(root))
            return
        entries = db_list_entries(root)
        typer.echo("PROVIDER\tSTATUS\tDBTYPE\tDESCRIPTION")
        for entry in entries:
            status = f"installed ({entry.records})" if entry.installed else "available"
            typer.echo(f"{entry.name}\t{status}\t{entry.dbtype}\t{entry.description}")

    dispatch(run)


def register_db_command(app: typer.Typer) -> None:
    """Attach the `db` command group (install, fetch, list) to the CLI app."""
    db_app = typer.Typer(
        help="Database acquisition and maintenance (provider fetch, verified local install).",
        no_args_is_help=True,
    )
    db_app.command("install")(db_install_command)
    db_app.command("build")(db_build_command)
    db_app.command("fetch")(db_fetch_command)
    db_app.command("list")(db_list_command)
    db_app.command("outdated")(db_outdated_command)
    db_app.command("search")(db_search_command)
    app.add_typer(db_app, name="db")
