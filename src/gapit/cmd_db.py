"""The `gapit db` command group: build, fetch, list, install.

- ``db build NAME FASTA``: build a custom gapit-native database from a
  user-supplied FASTA (+ optional TSV metadata) — the whole acquisition
  pipeline without a provider (:mod:`gapit.cmd_db_build`).
- ``db fetch [NAME]``: acquire provider database(s) under the datadir. With
  no NAME the default set (``DEFAULT_DBS``: card, vfdb) installs in order —
  each from its BUNDLED SNAPSHOT (Wave G) when one resolves, else over the
  network; ``--from-source`` forces the upstream download even when a
  snapshot exists. One JSON receipt line per database on stdout.
- ``db list``: show every known provider with its installed state.
- ``db outdated`` / ``db search``: read-only queries over the installed
  databases (:mod:`gapit.cmd_db_outdated`, :mod:`gapit.cmd_db_search`).
- ``db install``: checksum-verified LOCAL-FILE installation
  (:mod:`gapit.cmd_db_install`, split off in Wave G at the 250 LOC ceiling).
"""

from collections.abc import Iterator
from datetime import UTC, datetime
from pathlib import Path
from typing import Annotated, Literal

import typer
from pydantic import BaseModel, ConfigDict, Field

from gapit import config
from gapit.cmd_db_build import db_build_command
from gapit.cmd_db_install import db_install_command
from gapit.cmd_db_outdated import db_outdated_command
from gapit.cmd_db_search import db_search_command
from gapit.dispatch import Datadir, dispatch
from gapit.errors import UsageError
from gapit.providers import REGISTRY
from gapit.providers.common import Dbtype, fetch_provider
from gapit.records import read_manifest

# Bare `gapit db fetch` installs these, in order — the two providers whose
# snapshots ship inside the wheel (Wave G: zero-network bootstrap).
DEFAULT_DBS = ("card", "vfdb")


class ProviderReceipt(BaseModel, frozen=True):
    """One-line JSON success receipt for `db fetch` (no biological data)."""

    db: str
    records: int
    dbtype: Dbtype
    destination: str


class DbListEntry(BaseModel, frozen=True):
    """One provider row in the gapit.dblist/1 listing document."""

    name: str
    description: str
    dbtype: str
    installed: bool
    records: int | None = None


class DbListDocument(BaseModel, frozen=True):
    """gapit.dblist/1 — `gapit db list --json` output.

    A CLI listing, deliberately local to this module: NOT registered in
    `gapit schema` (the public schema surface stays unchanged this wave).
    """

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.dblist/1"] = Field(default="gapit.dblist/1", alias="schema")
    providers: tuple[DbListEntry, ...]


def perform_fetch(
    name: str | None,
    datadir: Path | None,
    *,
    force: bool = False,
    from_source: bool = False,
    quiet: bool = True,
    debug: bool = False,
) -> Iterator[ProviderReceipt]:
    """Fetch provider database(s) into <datadir>/NAME — the shared CLI + MCP
    path. NAME None installs every database in DEFAULT_DBS order; each
    receipt yields as its install completes (streaming, like the CLI's
    per-db stdout lines).
    """

    def fetch_one(provider_name: str) -> ProviderReceipt:
        provider = REGISTRY.get(provider_name)
        if provider is None:
            raise UsageError(
                f"unknown provider: {provider_name} (available: {', '.join(sorted(REGISTRY))})",
                code="USAGE_ERROR",
                context={"provider": provider_name},
            )
        db_dir = config.ensure_datadir(datadir) / provider_name
        manifest = fetch_provider(
            provider,
            db_dir,
            fetched_at=datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
            force=force,
            quiet=quiet,
            from_source=from_source,
            debug=debug,
        )
        return ProviderReceipt(
            db=provider_name,
            records=manifest.n_records,
            dbtype=manifest.dbtype,
            destination=str(db_dir),
        )

    for provider_name in (name,) if name is not None else DEFAULT_DBS:
        yield fetch_one(provider_name)


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


def db_list_entries(root: Path) -> list[DbListEntry]:
    """Provider rows for the gapit.dblist/1 listing — shared CLI + MCP path."""
    entries: list[DbListEntry] = []
    for provider_name in sorted(REGISTRY):
        provider = REGISTRY[provider_name]
        manifest_path = root / provider_name / "gapit-manifest.json"
        installed = manifest_path.is_file()
        entries.append(
            DbListEntry(
                name=provider_name,
                description=provider.description,
                dbtype=provider.dbtype,
                installed=installed,
                records=read_manifest(manifest_path).n_records if installed else None,
            )
        )
    return entries


def db_list_json(root: Path) -> str:
    """The gapit.dblist/1 document JSON — shared CLI + MCP path."""
    entries = db_list_entries(root)
    return DbListDocument(providers=tuple(entries)).model_dump_json(
        indent=2, by_alias=True, exclude_none=True
    )


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
