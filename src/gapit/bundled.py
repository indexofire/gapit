"""Bundled databases: license-clean content shipped inside the wheel.

A bundled database is a directory ``gapit/data/dbs/<name>/`` holding
``sequences``, a ``bundled.json`` metadata file, and kind-dependent
companions: gene bundles an optional ``typing.json`` and ``floors.json``
(gapit.floors/1), cluster bundles a REQUIRED ``features.json``
(gapit.features/1) plus their typing document. Materialization replays the
standard build pipelines into the user's datadir — the gene path through
:func:`gapit.db_build_ops.perform_build`, the cluster path through
:func:`gapit.bundled_cluster.materialize_cluster_bundle` (which rebuilds
from the committed locus FASTA + resolved feature table via
:func:`gapit.clusterbuild.build_cluster_database`) — deterministic, zero
network, so a bundled database is ready at install time with no
``db fetch``.

Only public-domain / permissive content may live here (the wheel is MIT):
PD, Apache-2.0, BSD, MIT — never GPL or non-commercial terms; everything
else stays download-on-fetch in the provider catalog. The three provider
snapshots (ncbi, resfinder, ecoh) are point-in-time copies
of exactly what ``db fetch`` would install (gapit/v1 headers, produced by
the provider pipeline itself); the typed cluster snapshot vp_ok is a pure
built-in (no provider behind it) carrying the resolved gapit.features/1
table a completed /2 build produced; ``db fetch <name>`` remains the
freshness path for the provider-backed bundled names and overwrites the
datadir copy with the latest upstream content. The gene
build pipeline re-headers any bundled FASTA into gapit/v1 like a custom
build, so both header styles round-trip.

Bundling is the exception, not the pattern, and every new bundle clears
two gates in order (the full policy lives in AGENTS.md §6): the license
gate above, then the size gate — the whole bundled root stays within a
30 MiB wheel budget, mechanically enforced by tests/test_bundled.py.
Past the budget, split a separate data wheel or keep the database
fetch-on-demand instead of inflating the ceiling.
"""

from pathlib import Path
from typing import Literal

from pydantic import BaseModel

from gapit import config
from gapit.bundled_cluster import materialize_cluster_bundle
from gapit.db_build_ops import perform_build
from gapit.db_info_models import DbGuidance
from gapit.errors import DatabaseError
from gapit.fasta import iter_fasta
from gapit.proctools import note

BUNDLED_ROOT = Path(__file__).resolve().parent / "data" / "dbs"


class BundledMetadata(BaseModel, frozen=True):
    """Parsed ``bundled.json`` — the listing-facing facts for one bundle.

    ``snapshotted`` is the ISO date the wheel copy was produced from
    upstream (provenance for humans; ``db fetch`` refreshes past it).
    ``kind`` routes materialization (gene vs cluster pipeline; cluster
    bundles additionally require ``features.json``). ``license`` and
    ``note`` carry content provenance into the materialized manifest
    (cluster bundles); they default to empty so the gene bundles'
    metadata stays minimal. ``guidance`` is what `gapit db info` reports:
    the literature references, the intended usage recipe, and the caveats
    a caller must not skip.
    """

    name: str
    description: str
    vendor: str
    dbtype: Literal["nucl", "prot"]
    kind: Literal["gene", "cluster"] = "gene"
    snapshotted: str
    license: str = ""
    note: str = ""
    guidance: DbGuidance = DbGuidance()


class BundledDatabase(BaseModel, frozen=True):
    """One bundled database directory plus its parsed metadata."""

    metadata: BundledMetadata
    path: Path

    @property
    def name(self) -> str:
        return self.metadata.name

    @property
    def sequences_path(self) -> Path:
        return self.path / "sequences"

    @property
    def typing_path(self) -> Path | None:
        """The bundled typing spec, when the database ships one."""
        typing = self.path / "typing.json"
        return typing if typing.is_file() else None

    @property
    def floors_path(self) -> Path | None:
        """The bundled gapit.floors/1 sidecar, when the database ships one."""
        floors = self.path / "floors.json"
        return floors if floors.is_file() else None

    @property
    def features_path(self) -> Path | None:
        """The bundled gapit.features/1 table — required (and checked at
        discovery) for cluster bundles, absent from gene bundles."""
        features = self.path / "features.json"
        return features if features.is_file() else None


def bundled_databases() -> list[BundledDatabase]:
    """Every bundled database, sorted by name. A directory qualifies when it
    carries both ``sequences`` and ``bundled.json`` — plus ``features.json``
    when the metadata declares ``kind: cluster`` (a cluster bundle without
    its feature table cannot materialize; typed error, no silent
    degradation). Malformed metadata propagates, and a wheel without the
    data dir lists nothing."""
    if not BUNDLED_ROOT.is_dir():
        return []
    databases: list[BundledDatabase] = []
    for child in sorted(BUNDLED_ROOT.iterdir()):
        metadata_path = child / "bundled.json"
        if not (child.is_dir() and (child / "sequences").is_file() and metadata_path.is_file()):
            continue
        try:
            metadata = BundledMetadata.model_validate_json(
                metadata_path.read_text(encoding="utf-8")
            )
        except (OSError, UnicodeDecodeError, ValueError) as exc:
            raise DatabaseError(
                f"bundled metadata is malformed: {metadata_path}: {exc}",
                code="BUNDLED_MALFORMED",
                context={"db": child.name, "file": str(metadata_path)},
            ) from exc
        if metadata.kind == "cluster" and not (child / "features.json").is_file():
            raise DatabaseError(
                f"bundled cluster database {metadata.name} lacks features.json",
                code="BUNDLED_MALFORMED",
                context={"db": metadata.name},
            )
        databases.append(BundledDatabase(metadata=metadata, path=child))
    return databases


def bundled_names() -> list[str]:
    """Sorted names of every bundled database."""
    return [database.name for database in bundled_databases()]


def find_bundled(name: str) -> BundledDatabase | None:
    """The bundled database called ``name``, or None when none ships."""
    return next((database for database in bundled_databases() if database.name == name), None)


def materialize_bundled(name: str, datadir: Path, *, quiet: bool = True) -> None:
    """Build a bundled database into ``<datadir>/name`` via its kind's
    standard pipeline (records.jsonl/sequences for gene bundles, locus
    FASTA + features.json for cluster bundles; BLAST index, typing copy,
    manifest stamped ``source: "bundled"`` — written last, as always).
    Zero network: every input travels inside the wheel.

    An already-materialized database (manifest present) is left untouched —
    callers short-circuit on the manifest, and each path's
    DB_ALREADY_EXISTS guard turns a racing double build into a typed error.
    """
    bundled = find_bundled(name)
    if bundled is None:
        raise DatabaseError(
            f"unknown bundled database: {name}",
            code="BUNDLED_NOT_FOUND",
            context={"db": name},
        )
    if bundled.metadata.kind == "cluster":
        materialize_cluster_bundle(bundled, datadir, quiet=quiet)
        return
    n_records = sum(1 for _ in iter_fasta(bundled.sequences_path))
    note(quiet, f"materializing bundled database {name} ({n_records} records) into {datadir}")
    perform_build(
        name,
        bundled.sequences_path,
        None,
        None,
        bundled.metadata.description,
        datadir,
        force=False,
        warn=lambda message: note(quiet, message),
        quiet=True,
        typing=bundled.typing_path,
        floors=bundled.floors_path,
        source="bundled",
    )


def resolve_screen_datadir(cli_value: Path | None, db_name: str) -> Path:
    """``resolve_datadir`` for the screening lookup — except that a MISSING
    datadir is bootstrapped when (and only when) the requested database is
    bundled: the zero-network first-use path must work on a fresh machine
    where ``~/.local/share/gapit/db`` does not exist yet. Every other
    missing-datadir case keeps the typed DATADIR_NOT_FOUND error."""
    try:
        return config.resolve_datadir(cli_value)
    except DatabaseError as exc:
        if exc.code != "DATADIR_NOT_FOUND" or find_bundled(db_name) is None:
            raise
        return config.ensure_datadir(cli_value)
