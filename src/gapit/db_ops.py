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
from gapit.bundled import BundledDatabase, bundled_databases
from gapit.db import Database, discover_databases, mol_type
from gapit.errors import UsageError
from gapit.fasta import iter_fasta
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

# PROVIDER cell for installed databases the registry does not know (db
# build / db install / abricate-style directories).
LOCAL_VENDOR = "local"


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
    # "local" marks a datadir-discovered database the registry does not know
    # (db build / db install product); "bundled" marks a database shipped in
    # the wheel (gapit.bundled) whether or not it is materialized yet.
    # Absent = registry catalog entry — the additive gapit.dblist/1 rule:
    # pre-existing entries stay byte-identical.
    source: Literal["local", "bundled"] | None = None


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


def _detected_dbtype(database: Database) -> str:
    """dbtype for a manifest-less database directory: BLAST index suffix
    first (``.nin`` nucl / ``.pin`` prot), else abricate's ``mol_type``
    heuristic over the sequence letters (the make_blast_db fallback)."""
    sequences = database.sequences_path
    for suffix, dbtype in ((".nin", "nucl"), (".pin", "prot")):
        if sequences.with_name(f"{sequences.name}{suffix}").exists():
            return dbtype
    letters = "".join(record.sequence for record in iter_fasta(sequences))
    return mol_type(letters)


def _extra_entries(root: Path, bundled: frozenset[str]) -> list[DbListEntry]:
    """Rows for installed databases neither the registry nor the wheel knows
    — `db build` products, checksum-installed or abricate-style directories.
    Bundled names are covered by their own rows above. Manifest-less
    directories count their FASTA records and detect dbtype from the index;
    DESCRIPTION comes from the manifest ``note`` when one is set. A malformed
    manifest propagates (the db.py discovery contract)."""
    entries: list[DbListEntry] = []
    for database in discover_databases(root):
        if database.name in REGISTRY or database.name in bundled:
            continue
        manifest_path = database.path / "gapit-manifest.json"
        manifest = read_manifest(manifest_path) if manifest_path.is_file() else None
        entries.append(
            DbListEntry(
                name=database.name,
                vendor=LOCAL_VENDOR,
                description=(manifest.note or "") if manifest else "",
                dbtype=manifest.dbtype if manifest else _detected_dbtype(database),
                installed=True,
                records=(
                    manifest.n_records
                    if manifest
                    else sum(1 for _ in iter_fasta(database.sequences_path))
                ),
                kind=database.kind,
                source="local",
            )
        )
    return entries


def _bundled_entries(bundled: Sequence[BundledDatabase], root: Path) -> list[DbListEntry]:
    """Rows for wheel-shipped databases: metadata (vendor, description,
    dbtype) from ``bundled.json``; installed state and record count from a
    materialized manifest under the datadir. A bundle whose directory exists
    but never finished materializing (no manifest) lists as not installed —
    the next screen/setupdb rebuilds it. A bundled name that is also a
    registry provider (ncbi, resfinder, ...) appears ONLY here: the bundled
    row replaces its registry row so one name is always one row."""
    entries: list[DbListEntry] = []
    for database in bundled:
        manifest_path = root / database.name / "gapit-manifest.json"
        manifest = read_manifest(manifest_path) if manifest_path.is_file() else None
        entries.append(
            DbListEntry(
                name=database.name,
                vendor=database.metadata.vendor,
                description=database.metadata.description,
                dbtype=manifest.dbtype if manifest else database.metadata.dbtype,
                installed=manifest is not None,
                records=manifest.n_records if manifest else None,
                kind=manifest.kind if manifest else "gene",
                source="bundled",
            )
        )
    return entries


def db_list_entries(root: Path) -> list[DbListEntry]:
    """Database rows for the gapit.dblist/1 listing — shared CLI + MCP path.

    Registry providers first (alphabetical, bundled names skipped — their
    rows come from the bundled section), then bundled databases by name,
    then datadir-discovered extras by name, so local ``db build``/``db
    install`` databases are never hidden.
    """
    bundled = tuple(bundled_databases())
    bundled_names = frozenset(database.name for database in bundled)
    entries: list[DbListEntry] = []
    for provider_name in sorted(REGISTRY):
        if provider_name in bundled_names:
            continue
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
    return entries + _bundled_entries(bundled, root) + _extra_entries(root, bundled_names)


def db_list_status(entry: DbListEntry) -> str:
    """STATUS cell text shared by the TSV and the rich table — one source so
    the two renderings can never diverge. Bundled-but-not-materialized rows
    read ``bundled`` (ready to materialize on first use); everything else
    uninstalled reads ``available`` (fetch-time download)."""
    if entry.installed:
        return f"installed ({entry.records})"
    return "bundled" if entry.source == "bundled" else "available"


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
