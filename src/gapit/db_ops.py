"""Database use-cases shared by the CLI and MCP: fetch and list.

``perform_fetch`` installs database(s) under the datadir (every provider
downloads from upstream at fetch time); the ``db_list_*`` callables build
the gapit.dblist/1 database listing. The typer commands
(:mod:`gapit.cmd_db`) and the MCP tools (:mod:`gapit.mcp_tools`) are thin
callers — this module owns the behavior and imports no CLI plumbing.
"""

from collections.abc import Iterator, Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field
from rich.table import Table
from rich.text import Text

from gapit import config
from gapit.errors import UsageError
from gapit.providers import REGISTRY
from gapit.providers.cluster_common import ClusterProvider, fetch_cluster_provider
from gapit.providers.common import Dbtype, fetch_provider
from gapit.records import read_manifest

# `gapit db fetch all` installs these, in order — the two headline
# databases (AMR + virulence). Like every provider they download on fetch:
# card (McMaster) and vfdb (CC BY-NC) licenses forbid redistribution.
DEFAULT_DBS = ("card", "vfdb")

# The literal NAME that selects DEFAULT_DBS; not a provider name (REGISTRY
# holds 19 names, none of which is "all").
FETCH_ALL = "all"


class ProviderReceipt(BaseModel, frozen=True):
    """One-line JSON success receipt for `db fetch` (no biological data)."""

    db: str
    records: int
    dbtype: Dbtype
    destination: str


class DbListEntry(BaseModel, frozen=True):
    """One database row in the gapit.dblist/1 listing document."""

    name: str
    # Upstream maintainer organisation; distinct from name (the --db value).
    vendor: str
    description: str
    dbtype: str
    installed: bool
    records: int | None = None
    # Installed providers inherit their manifest's kind (cluster for a
    # GBK/GFF-built db shadowing a provider name); the default is gene.
    kind: str = "gene"
    # Upstream content license when the provider pins one (card, vfdb,
    # ecoli_vf, kaptive); omitted otherwise.
    license: str | None = None


class DbListDocument(BaseModel, frozen=True):
    """gapit.dblist/1 — `gapit db list --json` output.

    A CLI listing, deliberately local to this module: NOT registered in
    `gapit schema` (the public schema surface stays unchanged this wave).
    """

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.dblist/1"] = Field(default="gapit.dblist/1", alias="schema")
    providers: tuple[DbListEntry, ...]


def perform_fetch(
    name: str,
    datadir: Path | None,
    *,
    force: bool = False,
    quiet: bool = True,
    debug: bool = False,
) -> Iterator[ProviderReceipt]:
    """Fetch database(s) into <datadir>/NAME — the shared CLI + MCP
    path. NAME ``all`` installs every database in DEFAULT_DBS order; each
    receipt yields as its install completes (streaming, like the CLI's
    per-db stdout lines).
    """

    def fetch_one(provider_name: str) -> ProviderReceipt:
        provider = REGISTRY.get(provider_name)
        if provider is None:
            raise UsageError(
                f"unknown database: {provider_name} (available: {', '.join(sorted(REGISTRY))})",
                code="USAGE_ERROR",
                context={"db": provider_name},
            )
        db_dir = config.ensure_datadir(datadir) / provider_name
        fetched_at = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
        if isinstance(provider, ClusterProvider):
            manifest = fetch_cluster_provider(
                provider,
                db_dir,
                fetched_at=fetched_at,
                force=force,
                quiet=quiet,
                debug=debug,
            )
        else:
            manifest = fetch_provider(
                provider,
                db_dir,
                fetched_at=fetched_at,
                force=force,
                quiet=quiet,
                debug=debug,
            )
        return ProviderReceipt(
            db=provider_name,
            records=manifest.n_records,
            dbtype=manifest.dbtype,
            destination=str(db_dir),
        )

    for provider_name in DEFAULT_DBS if name == FETCH_ALL else (name,):
        yield fetch_one(provider_name)


def db_list_entries(root: Path) -> list[DbListEntry]:
    """Database rows for the gapit.dblist/1 listing — shared CLI + MCP path."""
    entries: list[DbListEntry] = []
    for provider_name in sorted(REGISTRY):
        provider = REGISTRY[provider_name]
        manifest_path = root / provider_name / "gapit-manifest.json"
        installed = manifest_path.is_file()
        manifest = read_manifest(manifest_path) if installed else None
        entries.append(
            DbListEntry(
                name=provider_name,
                vendor=provider.vendor,
                description=provider.description,
                dbtype=provider.dbtype,
                installed=installed,
                records=manifest.n_records if manifest else None,
                kind=manifest.kind if manifest else provider.kind,
                license=provider.license,
            )
        )
    return entries


def db_list_status(entry: DbListEntry) -> str:
    """STATUS cell text shared by the TSV and the rich table — one source so
    the two renderings can never diverge."""
    return f"installed ({entry.records})" if entry.installed else "available"


def db_list_table(entries: Sequence[DbListEntry]) -> Table:
    """Rich rendering of the database listing for interactive terminals.

    Same rows and STATUS texts as the TSV (NAME, PROVIDER, STATUS, DBTYPE,
    DESCRIPTION); the command prints it only when stdout is a TTY, so pipes
    and redirects keep the byte-stable TSV. Description cells are literal
    ``Text`` so upstream prose can never parse as rich markup.
    """
    table = Table(title="Databases")
    table.add_column("Name", style="cyan")
    table.add_column("Provider")
    table.add_column("Status")
    table.add_column("DBTYPE", style="yellow")
    table.add_column("Description")
    for entry in entries:
        table.add_row(
            entry.name,
            entry.vendor,
            Text(db_list_status(entry), style="green" if entry.installed else "dim"),
            entry.dbtype,
            Text(entry.description),
        )
    return table


def db_list_json(root: Path) -> str:
    """The gapit.dblist/1 document JSON — shared CLI + MCP path."""
    entries = db_list_entries(root)
    return DbListDocument(providers=tuple(entries)).model_dump_json(
        indent=2, by_alias=True, exclude_none=True
    )
