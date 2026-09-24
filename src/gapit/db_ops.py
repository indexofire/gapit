"""Database provider use-cases shared by the CLI and MCP: fetch and list.

``perform_fetch`` installs provider database(s) under the datadir (bundled
snapshot first, ``from_source`` forces upstream); the ``db_list_*``
callables build the gapit.dblist/1 provider listing. The typer commands
(:mod:`gapit.cmd_db`) and the MCP tools (:mod:`gapit.mcp_tools`) are thin
callers — this module owns the behavior and imports no CLI plumbing.
"""

from collections.abc import Iterator
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from gapit import config
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
