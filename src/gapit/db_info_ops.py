"""The `db info` use-case: one database's provenance, literature, and caveats.

``db list`` answers "which databases exist and are they installed"; ``db
info`` answers "what is this database, where does it come from, how may I use
it, and where does it mislead". Everything it reports already exists in one of
four places, so the collector merges them instead of inventing a fifth:

- the **provider catalog** (:data:`gapit.providers.REGISTRY`) — description,
  vendor, license;
- the wheel's **``bundled.json``** (:mod:`gapit.bundled`) — the curated
  literature references, usage recipe, and caveats a bundled database ships
  with (``gapit.db_info_models.DbGuidance``);
- the installed **manifest** — what is actually on disk (record count, kind,
  dbtype, sha256, fetch time, typing schema, tool versions);
- an optional **``info.json``** sidecar in a database directory, for
  user-built databases that want to carry their own references/caveats.

The result is the ``gapit.dbinfo/1`` document (``db info --json``) or the
same facts as a stable ``FIELD<TAB>VALUE`` table (the default), list-valued
fields repeating their key — the agent-first shape ``db list`` established.
"""

from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from gapit import config
from gapit.bundled import BundledDatabase, find_bundled
from gapit.db_info_models import DbGuidance, DbReference
from gapit.db_ops import LOCAL_VENDOR, db_list_entries
from gapit.errors import DatabaseError
from gapit.formats.csvutil import tsv_cell
from gapit.providers import REGISTRY
from gapit.records import Manifest, read_manifest

INFO_FILENAME = "info.json"


class DbBuildInfo(BaseModel, frozen=True):
    """How the installed artifacts were produced (manifest provenance)."""

    tool: str
    fetched_at: str
    makeblastdb_version: str = ""
    minimap2_version: str = ""


class DbInfoDocument(BaseModel, frozen=True):
    """gapit.dbinfo/1 — one database's identity, provenance, and guidance.

    ``status`` mirrors ``db list``: ``installed`` (manifest on disk),
    ``bundled`` (ships in the wheel, materializes on first use), or
    ``available`` (fetch-time download). Provenance fields are omitted rather
    than null-filled when a source does not carry them, so a bundled panel and
    a fetched provider produce differently-shaped-but-honest documents.
    """

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.dbinfo/1"] = Field(default="gapit.dbinfo/1", alias="schema")
    name: str
    description: str
    vendor: str
    kind: Literal["gene", "cluster"]
    dbtype: Literal["nucl", "prot"]
    status: Literal["installed", "bundled", "available"]
    records: int | None = None
    destination: str | None = None
    license: str | None = None
    snapshot: str | None = None
    typing_schema: str | None = None
    source_urls: tuple[str, ...] = ()
    note: str | None = None
    build: DbBuildInfo | None = None
    references: tuple[DbReference, ...] = ()
    caveats: tuple[str, ...] = ()
    usage: tuple[str, ...] = ()


def _read_guidance(directory: Path) -> DbGuidance:
    """The optional ``info.json`` sidecar of one database directory.

    A user-built database can carry its own literature and caveats; a broken
    sidecar is a typed error rather than silently ignored provenance.
    """
    path = directory / INFO_FILENAME
    if not path.is_file():
        return DbGuidance()
    try:
        return DbGuidance.model_validate_json(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        raise DatabaseError(
            f"malformed database info sidecar: {path}: {exc}",
            code="DBINFO_MALFORMED",
            context={"file": str(path)},
        ) from exc


def _bundled_guidance(bundled: BundledDatabase | None) -> DbGuidance:
    """The wheel-declared literature/usage/caveats, or nothing."""
    if bundled is None:
        return DbGuidance()
    return bundled.metadata.guidance


def _default_usage(name: str, kind: str, typed: bool) -> tuple[str, ...]:
    """The workflow every database of this kind participates in.

    Declared guidance wins; this keeps every database's page runnable without
    each provider restating the same two commands.
    """
    if kind == "cluster":
        return (f"gapit screen -d {name} [FILE]... --format json",)
    if typed:
        return (
            f"gapit screen -d {name} [FILE]... -o RESULT.tsv",
            "gapit typing RESULT.tsv",
        )
    return (f"gapit screen -d {name} [FILE]... --format json",)


def _provider_guidance(name: str) -> DbGuidance:
    """Guidance declared by the provider module (empty for local databases)."""
    provider = REGISTRY.get(name)
    return provider.guidance if provider is not None else DbGuidance()


def _build_info(manifest: Manifest) -> DbBuildInfo:
    return DbBuildInfo(
        tool=f"{manifest.tool.name} {manifest.tool.version}",
        fetched_at=manifest.fetched_at,
        makeblastdb_version=manifest.makeblastdb_version,
        minimap2_version=manifest.minimap2_version,
    )


def collect_db_info(name: str, datadir: Path | None) -> DbInfoDocument:
    """Assemble one database's ``gapit.dbinfo/1`` document.

    Resolution order for the identity fields: the installed manifest (what is
    on disk) wins over the bundled metadata, which wins over the provider
    catalog; a name that is none of those is a ``DATABASE_NOT_FOUND``. The
    datadir bootstraps on a missing root (the ``db list`` rule) so ``db info``
    works on a fresh install before anything is fetched.
    """
    root = config.ensure_datadir(datadir)
    entry = next((row for row in db_list_entries(root) if row.name == name), None)
    if entry is None:
        raise DatabaseError(
            f"unknown database: {name} (see: gapit db list)",
            code="DATABASE_NOT_FOUND",
            context={"db": name},
        )
    directory = root / name
    bundled = find_bundled(name)
    manifest_path = directory / "gapit-manifest.json"
    manifest = read_manifest(manifest_path) if manifest_path.is_file() else None
    guidance = _bundled_guidance(bundled)
    if not guidance.references and not guidance.caveats:
        # Provider-declared guidance (the catalog's own literature notes).
        provider_guidance = _provider_guidance(name)
        guidance = provider_guidance
    if manifest is not None:
        sidecar = _read_guidance(directory)
        # The database's own sidecar wins over the wheel's declared guidance
        # (a user may override a bundled panel's caveats in their datadir).
        guidance = DbGuidance(
            references=sidecar.references or guidance.references,
            caveats=sidecar.caveats or guidance.caveats,
            usage=sidecar.usage or guidance.usage,
        )
    status: Literal["installed", "bundled", "available"]
    status = (
        "installed"
        if entry.installed
        else ("bundled" if entry.source == "bundled" else "available")
    )
    return DbInfoDocument(
        name=name,
        description=entry.description or (bundled.metadata.description if bundled else ""),
        vendor=entry.vendor
        if entry.vendor != LOCAL_VENDOR
        else (bundled.metadata.vendor if bundled else LOCAL_VENDOR),
        kind=manifest.kind if manifest else entry.kind,
        dbtype=manifest.dbtype if manifest else entry.dbtype,
        status=status,
        records=manifest.n_records if manifest else None,
        destination=str(directory) if manifest else None,
        # A materialized bundled panel keeps the wheel-declared licence (the
        # gene build does not carry one into the manifest).
        license=(manifest.license if manifest and manifest.license else entry.license),
        snapshot=bundled.metadata.snapshotted if bundled else None,
        typing_schema=manifest.typing_schema if manifest else None,
        source_urls=manifest.source_urls if manifest else (),
        note=(
            (manifest.note if manifest and manifest.note else None)
            or (bundled.metadata.note if bundled else None)
        ),
        build=_build_info(manifest) if manifest else None,
        references=guidance.references,
        caveats=guidance.caveats,
        usage=guidance.usage
        or _default_usage(
            name,
            manifest.kind if manifest else entry.kind,
            bool(manifest is not None and manifest.typing_schema),
        ),
    )


def render_db_info(document: DbInfoDocument) -> str:
    """The default ``FIELD<TAB>VALUE`` rendering, one row per fact.

    List-valued fields (``REFERENCE``, ``CAVEAT``, ``USAGE``, ``SOURCE_URL``)
    repeat their key, so the table stays two-column and stream-parseable while
    ordered lists survive verbatim; a reference's DOI/URL/note rows follow
    its own ``REFERENCE`` row. Empty fields are omitted — a bundled panel that
    is not yet materialized simply has no ``TOOL`` row.
    """
    rows: list[tuple[str, str]] = [
        ("NAME", document.name),
        ("DESCRIPTION", document.description),
        ("VENDOR", document.vendor),
        ("KIND", document.kind),
        ("DBTYPE", document.dbtype),
        ("STATUS", document.status),
        ("RECORDS", str(document.records) if document.records is not None else ""),
        ("DESTINATION", document.destination or ""),
        ("LICENSE", document.license or ""),
        ("SNAPSHOT", document.snapshot or ""),
        ("TYPING_SCHEMA", document.typing_schema or ""),
        ("NOTE", document.note or ""),
    ]
    if document.build is not None:
        rows.extend(
            [
                ("TOOL", document.build.tool),
                ("FETCHED_AT", document.build.fetched_at),
                ("MAKEBLASTDB_VERSION", document.build.makeblastdb_version),
                ("MINIMAP2_VERSION", document.build.minimap2_version),
            ]
        )
    rows.extend(("SOURCE_URL", url) for url in document.source_urls)
    # One block per reference, so a DOI/URL/note never floats free of the
    # citation it belongs to.
    for reference in document.references:
        rows.append(("REFERENCE", reference.citation))
        if reference.doi:
            rows.append(("REFERENCE_DOI", reference.doi))
        if reference.url:
            rows.append(("REFERENCE_URL", reference.url))
        if reference.note:
            rows.append(("REFERENCE_NOTE", reference.note))
    rows.extend(("CAVEAT", caveat) for caveat in document.caveats)
    rows.extend(("USAGE", step) for step in document.usage)
    return "".join(f"{field}\t{tsv_cell(value)}\n" for field, value in rows if value != "")
