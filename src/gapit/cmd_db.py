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
- ``db install``: checksum-verified LOCAL-FILE installation
  (:mod:`gapit.cmd_db_install`, split off in Wave G at the 250 LOC ceiling).
"""

from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path
from typing import Annotated, Literal

import typer
from pydantic import BaseModel, ConfigDict, Field

from gapit import config
from gapit.cmd_db_build import db_build_command
from gapit.cmd_db_install import db_install_command
from gapit.errors import DatabaseError, GapitError, UsageError, render_error
from gapit.providers import REGISTRY
from gapit.providers.common import Dbtype, fetch_provider
from gapit.records import read_manifest

# Bare `gapit db fetch` installs these, in order — the two providers whose
# snapshots ship inside the wheel (Wave G: zero-network bootstrap).
DEFAULT_DBS = ("card", "vfdb")

Datadir = Annotated[
    Path | None,
    typer.Option(
        "--datadir",
        help="Database directory (default: $GAPIT_DATADIR, then ~/.local/share/gapit/db).",
    ),
]


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


def _dispatch(action: Callable[[], None]) -> None:
    """Run a command body; failures render the gapit.error/1 envelope on
    stderr and exit with the documented code (mirrors cli._dispatch)."""
    try:
        action()
    except Exception as exc:
        typer.echo(render_error(exc), err=True)
        raise typer.Exit(code=exc.exit_code if isinstance(exc, GapitError) else 1) from exc


def _fetch_root(datadir: Path | None) -> Path:
    """Resolve the fetch datadir and mkdir it when absent (fresh-machine
    bootstrap, Wave E). The resolved path is recovered from resolve_datadir's
    DATADIR_NOT_FOUND context — config stays the single owner of resolution;
    only `db fetch` creates the root, every read path still demands it."""
    try:
        root = config.resolve_datadir(datadir)
    except DatabaseError as exc:
        if exc.code != "DATADIR_NOT_FOUND":
            raise
        root = Path(exc.context["datadir"])
    try:
        root.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise DatabaseError(
            f"cannot create datadir: {root}",
            code="DATADIR_CREATE_FAILED",
            context={"datadir": str(root)},
        ) from exc
    return root


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

    def fetch_one(provider_name: str) -> None:
        provider = REGISTRY.get(provider_name)
        if provider is None:
            raise UsageError(
                f"unknown provider: {provider_name} (available: {', '.join(sorted(REGISTRY))})",
                code="USAGE_ERROR",
                context={"provider": provider_name},
            )
        db_dir = _fetch_root(datadir) / provider_name
        manifest = fetch_provider(
            provider,
            db_dir,
            fetched_at=datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
            force=force,
            quiet=quiet,
            from_source=from_source,
            debug=debug,
        )
        typer.echo(
            ProviderReceipt(
                db=provider_name,
                records=manifest.n_records,
                dbtype=manifest.dbtype,
                destination=str(db_dir),
            ).model_dump_json()
        )

    def run() -> None:
        for provider_name in (name,) if name is not None else DEFAULT_DBS:
            fetch_one(provider_name)

    _dispatch(run)


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
        if as_json:
            document = DbListDocument(providers=tuple(entries))
            typer.echo(document.model_dump_json(indent=2, by_alias=True, exclude_none=True))
            return
        typer.echo("PROVIDER\tSTATUS\tDBTYPE\tDESCRIPTION")
        for entry in entries:
            status = f"installed ({entry.records})" if entry.installed else "available"
            typer.echo(f"{entry.name}\t{status}\t{entry.dbtype}\t{entry.description}")

    _dispatch(run)


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
    app.add_typer(db_app, name="db")
