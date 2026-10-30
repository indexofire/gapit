"""Database use-cases shared by the CLI and MCP: fetch and list.

``perform_fetch`` installs database(s) under the datadir (every provider
downloads from upstream at fetch time); the ``db_list_*`` callables build
the gapit.dblist/1 database listing. The CLI commands
(:mod:`gapit.cmd_db`) and the MCP tools (:mod:`gapit.mcp_tools`) are thin
callers — this module owns the behavior and imports no CLI plumbing.
"""

from collections.abc import Callable, Iterator, Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field
from rich import box
from rich.table import Table
from rich.text import Text

from gapit import config
from gapit.bundled import BundledDatabase, bundled_databases
from gapit.db import Database, discover_databases, mol_type
from gapit.db_info_models import DbGuidance
from gapit.errors import DatabaseError, UsageError
from gapit.fasta import iter_fasta
from gapit.progress import screen_progress
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
    dbtype: Literal["nucl", "prot"]
    installed: bool
    indexed: bool = False
    records: int | None = None
    # Installed providers inherit their manifest's kind (cluster for a
    # GBK/GFF-built db shadowing a provider name); the default is gene.
    kind: Literal["gene", "cluster"] = "gene"
    # Whether the database carries a typing scheme (typing.json). Gene-kind
    # databases designate through `gapit typing`; cluster-kind databases
    # designate inside `gapit screen` (the command refuses cluster dbs).
    typing: bool = False
    # Upstream content license when the provider pins one (card, vfdb,
    # ecoli_vf, kaptive); omitted otherwise.
    license: str | None = None
    # "local" marks a datadir-discovered database the registry does not know
    # (db build / db install product); "bundled" marks a database shipped in
    # the wheel (gapit.bundled) whether or not it is materialized yet.
    # Absent = registry catalog entry — the additive gapit.dblist/1 rule:
    # pre-existing entries stay byte-identical.
    source: Literal["local", "bundled"] | None = None


class DatadirDocument(BaseModel, frozen=True):
    """gapit.datadir/1 — the resolved datadir plus how it was chosen.

    A CLI introspection document (like :class:`DbListDocument`, local to this
    module): agents and container entrypoints can discover the effective path
    without guessing environment precedence.
    """

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.datadir/1"] = Field(default="gapit.datadir/1", alias="schema")
    datadir: str
    exists: bool
    source: Literal["cli", "env", "xdg", "default"]


class DbListDocument(BaseModel, frozen=True):
    """gapit.dblist/1 — `gapit db list --json` output.

    Registered as ``gapit schema dblist``; the schema is additive so agents can
    introspect the database catalog without a live datadir.
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
    emit: Callable[[ProviderReceipt], None] | None = None,
) -> Iterator[ProviderReceipt]:
    """Fetch database(s) into <datadir>/NAME — the shared CLI + MCP
    path. NAME ``all`` installs every database in DEFAULT_DBS order; each
    receipt yields as its install completes (streaming, like the CLI's
    per-db stdout lines). On an interactive stderr the batch rides the
    shared progress bar; under --quiet / non-TTY stderr the no-op shim
    keeps today's byte-identical behavior.
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

    names = DEFAULT_DBS if name == FETCH_ALL else (name,)
    with screen_progress(names, quiet=quiet, description="Fetching") as (bar, task_id):
        for provider_name in names:
            if task_id is not None:
                bar.describe(provider_name)
            receipt = fetch_one(provider_name)
            bar.advance()
            if emit is not None:
                with bar.suspend():
                    emit(receipt)
            yield receipt


def datadir_json(cli_value: Path | None) -> str:
    """The gapit.datadir/1 document for the resolved datadir (CLI + agents)."""
    path, source = config.datadir_source(cli_value)
    return DatadirDocument(datadir=str(path), exists=path.is_dir(), source=source).model_dump_json(
        by_alias=True
    )


def _detected_dbtype(database: Database) -> Literal["nucl", "prot"]:
    """dbtype for a manifest-less database directory: BLAST index suffix
    first (``.nin`` nucl / ``.pin`` prot), else abricate's ``mol_type``
    heuristic over the sequence letters (the make_blast_db fallback)."""
    sequences = database.sequences_path
    index_types: tuple[tuple[str, Literal["nucl", "prot"]], ...] = (
        (".nin", "nucl"),
        (".pin", "prot"),
    )
    for suffix, detected in index_types:
        if sequences.with_name(f"{sequences.name}{suffix}").exists():
            return detected
    letters = "".join(record.sequence for record in iter_fasta(sequences))
    return mol_type(letters)


def _local_listing_metadata(directory: Path) -> tuple[str, str]:
    """Description/vendor from a local ``info.json`` sidecar, when present.

    Local databases built with ``db build`` can carry the same ``info.json``
    guidance document used by ``db info``; its optional ``description`` and
    ``vendor`` fields feed the listing so curated panels are not anonymous
    ``local`` rows. A malformed sidecar is a typed error, never silently
    ignored provenance.
    """
    path = directory / "info.json"
    if not path.is_file():
        return "", ""
    try:
        guidance = DbGuidance.model_validate_json(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        raise DatabaseError(
            f"malformed database info sidecar: {path}: {exc}",
            code="DBINFO_MALFORMED",
            context={"file": str(path)},
        ) from exc
    return guidance.description, guidance.vendor


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
        info_description, info_vendor = _local_listing_metadata(database.path)
        entries.append(
            DbListEntry(
                name=database.name,
                vendor=info_vendor or LOCAL_VENDOR,
                description=info_description or ((manifest.note or "") if manifest else ""),
                dbtype=manifest.dbtype if manifest else _detected_dbtype(database),
                installed=True,
                indexed=_has_blast_index(database.sequences_path),
                records=(
                    manifest.n_records
                    if manifest
                    else sum(1 for _ in iter_fasta(database.sequences_path))
                ),
                kind=database.kind,
                typing=_has_typing(database.path),
                source="local",
            )
        )
    return entries


def _has_blast_index(sequences_path: Path) -> bool:
    """Whether either BLAST index family exists for this sequences file."""
    return any(
        (sequences_path.parent / f"{sequences_path.name}{suffix}").exists()
        for suffix in (".nin", ".pin")
    )


def _has_typing(path: Path) -> bool:
    """Whether a database directory carries a typing scheme."""
    return (path / "typing.json").is_file()


def _provider_has_typing(provider: object) -> bool:
    """Whether an unfetched provider would install a typing scheme (only
    cluster providers can declare the asset; gene providers never do)."""
    return isinstance(provider, ClusterProvider) and provider.typing is not None


def _bundled_entries(bundled: Sequence[BundledDatabase], root: Path) -> list[DbListEntry]:
    """Rows for wheel-shipped databases: metadata (vendor, description,
    dbtype) from ``bundled.json``; installed state and record count from a
    materialized manifest under the datadir. A bundle whose directory exists
    but never finished materializing (no manifest) lists as not installed —
    the next screen or `db setup` rebuilds it. A bundled name that is also a
    registry provider (ncbi, resfinder, ...) appears ONLY here: the bundled
    row replaces its registry row so one name is always one row."""
    entries: list[DbListEntry] = []
    for database in bundled:
        manifest_path = root / database.name / "gapit-manifest.json"
        installed = manifest_path.is_file() and (root / database.name / "sequences").is_file()
        manifest = read_manifest(manifest_path) if installed else None
        entries.append(
            DbListEntry(
                name=database.name,
                vendor=database.metadata.vendor,
                description=database.metadata.description,
                dbtype=manifest.dbtype if manifest else database.metadata.dbtype,
                installed=manifest is not None,
                indexed=(
                    _has_blast_index(root / database.name / "sequences")
                    if manifest is not None
                    else False
                ),
                records=manifest.n_records if manifest else None,
                kind=manifest.kind if manifest else database.metadata.kind,
                typing=(
                    _has_typing(root / database.name)
                    if manifest is not None
                    else database.typing_path is not None
                ),
                license=database.metadata.license or None,
                source="bundled",
            )
        )
    return entries


def db_list_entries(root: Path) -> list[DbListEntry]:
    """Database rows for the gapit.dblist/1 listing — shared CLI + MCP path.

    Rows are sorted globally by NAME so CLI TSV, the rich table, and the
    ``gapit.dblist/1`` JSON all share one deterministic order. Registry
    providers, bundled databases, and datadir-discovered local databases
    stay distinguishable by their additive ``source`` field and STATUS, not
    by position. Bundled names are skipped in the registry loop because
    their rows come from the bundled section; local ``db build``/``db
    install`` databases are appended and never hidden.
    """
    bundled = tuple(bundled_databases())
    bundled_names = frozenset(database.name for database in bundled)
    entries: list[DbListEntry] = []
    for provider_name in sorted(REGISTRY):
        if provider_name in bundled_names:
            continue
        provider = REGISTRY[provider_name]
        manifest_path = root / provider_name / "gapit-manifest.json"
        installed = manifest_path.is_file() and (root / provider_name / "sequences").is_file()
        manifest = read_manifest(manifest_path) if installed else None
        entries.append(
            DbListEntry(
                name=provider_name,
                vendor=provider.vendor,
                description=provider.description,
                dbtype=provider.dbtype,
                installed=installed,
                indexed=(
                    _has_blast_index(root / provider_name / "sequences") if installed else False
                ),
                records=manifest.n_records if manifest else None,
                kind=manifest.kind if manifest else provider.kind,
                typing=(
                    _has_typing(root / provider_name)
                    if installed
                    else _provider_has_typing(provider)
                ),
                license=provider.license,
            )
        )
    entries.extend(_bundled_entries(bundled, root))
    entries.extend(_extra_entries(root, bundled_names))
    return sorted(entries, key=lambda entry: entry.name)


def db_list_status(entry: DbListEntry) -> str:
    """STATUS cell text shared by the TSV and the rich table — one source so
    the two renderings can never diverge. Bundled-but-not-materialized rows
    read ``bundled`` (ready to materialize on first use); everything else
    uninstalled reads ``available`` (fetch-time download)."""
    if entry.installed:
        return f"installed ({entry.records})"
    return "bundled" if entry.source == "bundled" else "available"


def db_list_typing(entry: DbListEntry) -> str:
    """TYPING cell text: how the database designates, when it can.

    ``command`` — a gene-kind database: screen to a table, then run
    ``gapit typing RESULT.tsv``. ``integrated`` — a typed cluster database:
    the phenotype rides ``gapit screen`` itself (``gapit typing`` refuses
    cluster dbs). ``-`` — no typing.json, screen-only."""
    if not entry.typing:
        return "-"
    return "command" if entry.kind == "gene" else "integrated"


def db_list_table(entries: Sequence[DbListEntry]) -> Table:
    """Rich rendering of the database listing for interactive terminals.

    Same rows and STATUS/TYPING texts as the TSV (NAME, DESCRIPTION, STATUS,
    DBTYPE, TYPING, PROVIDER); the command prints it only when stdout is a
    TTY, so pipes and redirects keep the byte-stable TSV. Description cells
    are literal ``Text`` so upstream prose can never parse as rich markup.
    The box is SQUARE: every horizontal rule the same thin ``─`` weight (the
    rich default draws the top and head rules heavy).
    """
    table = Table(title="Databases", box=box.SQUARE)
    table.add_column("Name", style="cyan", no_wrap=True)
    table.add_column("Description")
    table.add_column("Status", no_wrap=True)
    table.add_column("DBTYPE", style="yellow", no_wrap=True)
    table.add_column("TYPING", style="magenta", no_wrap=True)
    table.add_column("Provider")
    for entry in entries:
        table.add_row(
            entry.name,
            Text(entry.description),
            Text(db_list_status(entry), style="green" if entry.installed else "dim"),
            entry.dbtype,
            db_list_typing(entry),
            entry.vendor,
        )
    return table


def db_list_json(root: Path) -> str:
    """The gapit.dblist/1 document JSON — shared CLI + MCP path."""
    entries = db_list_entries(root)
    return DbListDocument(providers=tuple(entries)).model_dump_json(
        indent=2, by_alias=True, exclude_none=True
    )
