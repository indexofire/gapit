"""Rich cluster locus-map models (prototype surface).

These are additive, optional report payloads: the existing TSV/MD and the
compact cluster/1-2 JSON stay unchanged unless a caller explicitly builds a
locus map.  The map is the shared input for quality decisions and external
rendering (agents/tools own the presentation layer).
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel

from gapit.cluster_alignment import AlignmentPiece, Deletion, Insertion
from gapit.gene_families import ContainmentEdge


class LocusMapGene(BaseModel, frozen=True):
    """One reference gene with screen evidence and family assignment."""

    gene_id: str
    start: int
    end: int
    strand: Literal["+", "-"]
    product: str
    family: str | None
    coverage_pct: float
    identity_pct: float
    verdict: Literal["present", "partial", "absent"]


class ExtraPanelGene(BaseModel, frozen=True):
    """A gene from another panel locus found inside the best locus region."""

    gene_id: str
    family: str
    source_locus: str
    contig: str
    qstart: int
    qend: int


class LocusQuality(BaseModel, frozen=True):
    """Prototype locus-quality decision (annotate-only; never gates typing)."""

    status: Literal["pass", "fail"]
    reasons: tuple[str, ...]
    min_coverage: float
    min_identity: float
    max_missing_genes: int
    max_extra_panel_genes: int
    coverage_pct: float
    identity_pct: float
    missing_genes: int
    partial_genes: int
    extra_panel_genes: int
    extra_bp: int
    containment_conflict: bool


class LocusMap(BaseModel, frozen=True):
    """Detailed evidence for one (best) locus."""

    locus: str
    label: str
    type: str
    length: int
    coverage_pct: float
    identity_pct: float
    genes: tuple[LocusMapGene, ...]
    pieces: tuple[AlignmentPiece, ...]
    insertions: tuple[Insertion, ...]
    deletions: tuple[Deletion, ...]
    missing: tuple[str, ...]
    partial: tuple[str, ...]
    extra_panel_genes: tuple[ExtraPanelGene, ...]
    containment: tuple[ContainmentEdge, ...]
    quality: LocusQuality
