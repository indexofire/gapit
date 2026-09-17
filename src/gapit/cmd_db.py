"""The `gapit db` command group: install, fetch, list.

- ``db install``: checksum-verified LOCAL-FILE installation. SOURCE must be a
  path to an existing regular file; no network, no provider IDs, no archives,
  no manifests — streaming SHA256 + atomic copy, kept deliberately independent
  of db.py and the screening pipeline.
- ``db fetch NAME``: acquire a named provider database over the network
  (Wave B provider registry) and build it under the datadir.
- ``db list``: show every known provider with its installed state.
"""

import hashlib
import os
import re
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Annotated, Literal

import typer
from pydantic import BaseModel, ConfigDict, Field

from gapit import config
from gapit.errors import DatabaseError, GapitError, InputError, UsageError, render_error
from gapit.providers import REGISTRY
from gapit.providers.common import Dbtype, fetch_provider
from gapit.records import read_manifest

_SHA256_SHAPE = re.compile(r"[0-9a-fA-F]{64}")
_CHUNK_BYTES = 1 << 20

Datadir = Annotated[
    Path | None,
    typer.Option(
        "--datadir",
        help="Database directory (default: $GAPIT_DATADIR, then ~/.local/share/gapit/db).",
    ),
]


class FetchReceipt(BaseModel, frozen=True):
    """One-line JSON success receipt printed to stdout (no biological data)."""

    destination: str
    sha256: str


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


def _parse_sha256(raw: str) -> str:
    """Lowercased digest if it is exactly 64 hex characters, else UsageError."""
    if not _SHA256_SHAPE.fullmatch(raw):
        raise UsageError(
            "--sha256 must be exactly 64 hex characters",
            code="USAGE_ERROR",
            context={"sha256": raw},
        )
    return raw.lower()


def install_verified(source: Path, target: Path, expected: str) -> FetchReceipt:
    """Copy source to target with a streaming SHA256 check; atomic on success.

    Bytes land in a NamedTemporaryFile in the target's parent (same
    filesystem) and are os.replace()d over the target only after the digest
    verifies, so any failure leaves an existing target untouched. The temp
    file is removed on every error path.
    """
    if not source.is_file():
        raise InputError(
            f"source file not found or unreadable: {source}",
            code="INPUT_NOT_FOUND",
            context={"source": str(source)},
        )
    temp_path: Path | None = None
    try:
        with source.open("rb") as src, NamedTemporaryFile(dir=target.parent, delete=False) as temp:
            temp_path = Path(temp.name)
            hasher = hashlib.sha256()
            while chunk := src.read(_CHUNK_BYTES):
                hasher.update(chunk)
                temp.write(chunk)
            actual = hasher.hexdigest()
            if actual != expected:
                raise InputError(
                    f"SHA256 mismatch for {source}",
                    code="CHECKSUM_MISMATCH",
                    context={
                        "source": str(source),
                        "expected": expected,
                        "actual": actual,
                    },
                )
        assert temp_path is not None
        os.replace(temp_path, target)
    except OSError as exc:
        raise InputError(
            f"cannot install {source} to {target}: {exc}",
            code="FILE_IO_ERROR",
            context={"source": str(source), "target": str(target)},
        ) from exc
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
    return FetchReceipt(destination=str(target), sha256=expected)


def db_install_command(
    source: Annotated[
        Path,
        typer.Argument(
            help="Local source file path (plain filesystem only; no URLs, no provider IDs)."
        ),
    ],
    sha256: Annotated[
        str,
        typer.Option("--sha256", help="Expected SHA256 digest of the source (64 hex chars)."),
    ],
    output: Annotated[
        Path,
        typer.Option(
            "--output", help="Destination path; replaced atomically only after verification."
        ),
    ],
) -> None:
    """Install a local file to --output after verifying its SHA256.

    VERIFIED LOCAL-FILE INSTALLATION ONLY: this never fetches over the
    network and knows nothing about database providers or sequence content —
    it checksum-verifies and atomically installs bytes. On success a one-line
    JSON receipt (destination, verified digest) is printed to stdout.
    """

    def run() -> None:
        expected = _parse_sha256(sha256)
        receipt = install_verified(source, output, expected)
        typer.echo(receipt.model_dump_json())

    _dispatch(run)


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
        str,
        typer.Argument(help="Provider name (see: gapit db list)."),
    ],
    datadir: Datadir = None,
    force: Annotated[
        bool,
        typer.Option("--force", help="Overwrite the database if it already exists."),
    ] = False,
    quiet: Annotated[bool, typer.Option("--quiet", help="Silence stderr diagnostics.")] = False,
) -> None:
    """Fetch and build a provider database into <datadir>/NAME (network)."""

    def run() -> None:
        provider = REGISTRY.get(name)
        if provider is None:
            raise UsageError(
                f"unknown provider: {name} (available: {', '.join(sorted(REGISTRY))})",
                code="USAGE_ERROR",
                context={"provider": name},
            )
        db_dir = _fetch_root(datadir) / name
        manifest = fetch_provider(
            provider,
            db_dir,
            fetched_at=datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
            force=force,
            quiet=quiet,
        )
        typer.echo(
            ProviderReceipt(
                db=name,
                records=manifest.n_records,
                dbtype=manifest.dbtype,
                destination=str(db_dir),
            ).model_dump_json()
        )

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
    db_app.command("fetch")(db_fetch_command)
    db_app.command("list")(db_list_command)
    app.add_typer(db_app, name="db")
