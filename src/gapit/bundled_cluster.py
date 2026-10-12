"""Kind-aware materialization for CLUSTER bundles (the bundled.py sibling).

A cluster bundle ships its build INPUTS committed: the locus FASTA plus the
``gapit.features/1`` table already carrying RESOLVED gene ids (the snapshot
was produced by a completed gapit.typing/2 build, so its rules reference
the resolved locus_tag-style ids). Materialization reconstructs
:class:`LocusSequence` objects from those artifacts and replays
:func:`gapit.clusterbuild.build_cluster_database` — the entry taking
pre-parsed loci, the SAME code a provider fetch runs — into the user's
datadir, deterministic and zero network.

Gene-id resolution deliberately does NOT re-run: the bundled feature table
is already resolved, and re-running would be a no-op at best. The typing
document IS validated against the reconstructed loci (mirroring the build
path's pre-write gate), and the screen-time /2 guard re-checks afterwards
— belt and braces, like the fetch path.
"""

from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING

from pydantic import ValidationError

from gapit.clusterbuild import build_cluster_database
from gapit.errors import DatabaseError
from gapit.fasta import iter_fasta
from gapit.gbfeatures import FeaturesDocument, LocusSequence
from gapit.proctools import note
from gapit.typing_models import (
    read_typing_document,
    typing_schema_of,
    validate_references,
)

if TYPE_CHECKING:  # bundled.py imports this module at runtime; no cycle.
    from gapit.bundled import BundledDatabase


def _malformed(db_name: str, reason: str) -> DatabaseError:
    """The broken-wheel-cluster-bundle failure: located, stable code."""
    return DatabaseError(
        f"bundled cluster database {db_name} is malformed ({reason})",
        code="BUNDLED_MALFORMED",
        context={"db": db_name, "reason": reason},
    )


def _reconstruct_loci(bundled: "BundledDatabase") -> list[LocusSequence]:
    """The bundle's committed artifacts back into LocusSequence objects.

    features.json supplies the feature tables (loci in document order), the
    FASTA the sequences (one bare-id record per locus); every locus must
    pair with the same-order same-id record — the reconstruction mirror of
    the build self-check, BUNDLED_MALFORMED on any drift."""
    features_path = bundled.path / "features.json"
    try:
        document = FeaturesDocument.model_validate_json(features_path.read_text(encoding="utf-8"))
    except ValidationError as exc:
        raise _malformed(bundled.name, f"unparseable features.json: {exc}") from exc
    records = list(iter_fasta(bundled.sequences_path))
    if len(records) != len(document.loci):
        raise _malformed(
            bundled.name,
            f"{len(document.loci)} loci but {len(records)} sequence records",
        )
    loci: list[LocusSequence] = []
    for locus, record in zip(document.loci, records, strict=True):
        if locus.id != record.id:
            raise _malformed(bundled.name, f"locus {locus.id!r} paired with record {record.id!r}")
        loci.append(LocusSequence(features=locus, sequence=record.sequence))
    return loci


def materialize_cluster_bundle(
    bundled: "BundledDatabase", datadir: Path, *, quiet: bool = True
) -> None:
    """Build a cluster-kind bundle into ``<datadir>/name`` through the
    cluster pipeline (locus FASTA, features.json, typing copy, BLAST index,
    manifest stamped ``source: "bundled"`` and written last, as always).
    The typing spec validates BEFORE any artifact is written, so a spec
    drifting from the committed table fails with no half-built database;
    an existing manifest is a typed DB_ALREADY_EXISTS (the racing
    double-build guard perform_build gives the gene path)."""
    name = bundled.name
    db_dir = datadir / name
    if (db_dir / "gapit-manifest.json").is_file():
        raise DatabaseError(
            f"won't overwrite existing database {name} (use --force)",
            code="DB_ALREADY_EXISTS",
            context={"db": name},
        )
    n_loci = sum(1 for _ in iter_fasta(bundled.sequences_path))
    note(quiet, f"materializing bundled database {name} ({n_loci} loci) into {datadir}")
    loci = _reconstruct_loci(bundled)
    typing_schema = ""
    if bundled.typing_path is not None:
        document = read_typing_document(bundled.typing_path)
        validate_references(
            document, FeaturesDocument(loci=tuple(locus.features for locus in loci))
        )
        typing_schema = typing_schema_of(bundled.typing_path)
    db_dir.mkdir(parents=True, exist_ok=True)
    build_cluster_database(
        db_dir,
        name=name,
        loci=loci,
        typing=bundled.typing_path,
        source_urls=("local",),
        fetched_at=datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        typing_schema=typing_schema,
        quiet=True,  # the materializing note above is the only stderr line (gene-path parity)
        content_license=bundled.metadata.license,
        content_note=bundled.metadata.note,
        source="bundled",
    )
