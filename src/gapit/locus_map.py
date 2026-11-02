"""Assemble a rich LocusMap from the compact screen report.

This is the prototype bridge between the existing LocusCall/typing output and
the detailed alignment evidence.  It deliberately does not change the compact
rendering: callers opt in by building the map, and quality is annotate-only.
"""

from __future__ import annotations

from collections.abc import Mapping

from gapit.cluster import GeneCall, LocusCall
from gapit.cluster_alignment import LocusAlignment, QuerySpan
from gapit.gbfeatures import FeaturesDocument, LocusFeatures
from gapit.gene_families import (
    GeneFamilyDocument,
    containment_for,
    gene_family_map,
    locus_family_sets,
)
from gapit.locus_map_models import (
    ExtraPanelGene,
    LocusMap,
    LocusMapGene,
    LocusQuality,
)


def _features_for(features: FeaturesDocument, locus: str) -> LocusFeatures:
    for entry in features.loci:
        if entry.id == locus:
            return entry
    raise ValueError(f"locus {locus!r} not found in features")


def _gene_calls(best: LocusCall) -> dict[str, GeneCall]:
    return {gene.gene_id: gene for gene in best.genes}


def _containment_conflict(
    best: LocusCall,
    families: GeneFamilyDocument | None,
    locus_calls: dict[str, LocusCall],
    *,
    min_coverage: float,
    min_identity: float,
) -> bool:
    """True when a superset locus is also supported and its diagnostic
    families have sample evidence."""
    if families is None:
        return False
    family_of = gene_family_map(families)
    for edge in containment_for(families, best.locus):
        if edge.subset != best.locus:
            continue
        superset = locus_calls.get(edge.superset)
        if superset is None:
            continue
        if superset.coverage_pct < min_coverage or superset.identity_pct < min_identity:
            continue
        diagnostic = set(edge.diagnostic)
        if any(
            family_of.get(gene.gene_id) in diagnostic and gene.verdict != "absent"
            for gene in superset.genes
        ):
            return True
    return False


def _overlap(left: tuple[str, int, int], right: QuerySpan) -> bool:
    return left[0] == right.contig and left[1] < right.qend and left[2] > right.qstart


def _extra_panel_genes(
    best: LocusCall,
    alignment: LocusAlignment,
    loci: tuple[LocusCall, ...],
    families: GeneFamilyDocument | None,
    gene_spans: Mapping[str, Mapping[str, tuple[QuerySpan, ...]]] | None,
) -> tuple[ExtraPanelGene, ...]:
    """Panel genes from other loci whose detected query span overlaps the
    best locus, excluding families already expected in the best locus."""
    if families is None or gene_spans is None:
        return ()
    family_of = gene_family_map(families)
    best_families = locus_family_sets(families).get(best.locus, frozenset())
    best_spans = [(piece.contig, piece.qstart, piece.qend) for piece in alignment.pieces]
    found: list[ExtraPanelGene] = []
    seen: set[tuple[str, str, str]] = set()
    for other in loci:
        if other.locus == best.locus:
            continue
        spans = gene_spans.get(other.locus, {})
        for gene in other.genes:
            if gene.verdict == "absent":
                continue
            family = family_of.get(gene.gene_id)
            if family is None or family in best_families:
                continue
            for span in spans.get(gene.gene_id, ()):
                if not any(_overlap(best_span, span) for best_span in best_spans):
                    continue
                key = (gene.gene_id, family, other.locus)
                if key in seen:
                    break
                seen.add(key)
                found.append(
                    ExtraPanelGene(
                        gene_id=gene.gene_id,
                        family=family,
                        source_locus=other.locus,
                        contig=span.contig,
                        qstart=span.qstart,
                        qend=span.qend,
                    )
                )
                break
    return tuple(sorted(found, key=lambda item: (item.source_locus, item.gene_id)))


def build_locus_map(
    best: LocusCall,
    alignment: LocusAlignment,
    features: FeaturesDocument,
    loci: tuple[LocusCall, ...],
    families: GeneFamilyDocument | None = None,
    gene_spans: Mapping[str, Mapping[str, tuple[QuerySpan, ...]]] | None = None,
    *,
    min_coverage: float = 80.0,
    min_identity: float = 90.0,
    max_missing_genes: int = 3,
    max_extra_panel_genes: int = 2,
) -> LocusMap:
    """Build the detailed map for one best locus (prototype thresholds)."""
    locus_features = _features_for(features, best.locus)
    calls = _gene_calls(best)
    family_of = gene_family_map(families) if families is not None else {}
    genes = tuple(
        LocusMapGene(
            gene_id=gene.gene_id,
            start=gene.start,
            end=gene.end,
            strand=gene.strand,
            product=gene.product,
            family=family_of.get(gene.gene_id),
            coverage_pct=calls[gene.gene_id].coverage_pct,
            identity_pct=calls[gene.gene_id].identity_pct,
            verdict=calls[gene.gene_id].verdict,
        )
        for gene in locus_features.genes
    )
    locus_calls = {call.locus: call for call in loci}
    conflict = _containment_conflict(
        best,
        families,
        locus_calls,
        min_coverage=min_coverage,
        min_identity=min_identity,
    )
    insertions = alignment.insertions
    extra_panel = _extra_panel_genes(best, alignment, loci, families, gene_spans)
    extra_bp = sum(item.length for item in insertions)
    missing = tuple(gene.gene_id for gene in genes if gene.verdict == "absent")
    partial = tuple(gene.gene_id for gene in genes if gene.verdict == "partial")
    reasons: list[str] = []
    if best.coverage_pct < min_coverage:
        reasons.append(f"coverage {best.coverage_pct:.2f} < {min_coverage:.2f}")
    if best.identity_pct < min_identity:
        reasons.append(f"identity {best.identity_pct:.2f} < {min_identity:.2f}")
    if len(missing) > max_missing_genes:
        reasons.append(f"missing genes {len(missing)} > {max_missing_genes}")
    if len(extra_panel) > max_extra_panel_genes:
        reasons.append(f"extra panel genes {len(extra_panel)} > {max_extra_panel_genes}")
    if conflict:
        reasons.append("containment conflict: subset locus with supported superset")
    quality = LocusQuality(
        status="pass" if not reasons else "fail",
        reasons=tuple(reasons),
        min_coverage=min_coverage,
        min_identity=min_identity,
        max_missing_genes=max_missing_genes,
        max_extra_panel_genes=max_extra_panel_genes,
        coverage_pct=best.coverage_pct,
        identity_pct=best.identity_pct,
        missing_genes=len(missing),
        partial_genes=len(partial),
        extra_panel_genes=len(extra_panel),
        extra_bp=extra_bp,
        containment_conflict=conflict,
    )
    return LocusMap(
        locus=best.locus,
        label=best.label,
        type=best.type,
        length=alignment.length,
        coverage_pct=best.coverage_pct,
        identity_pct=best.identity_pct,
        genes=genes,
        pieces=alignment.pieces,
        insertions=alignment.insertions,
        deletions=alignment.deletions,
        missing=missing,
        partial=partial,
        extra_panel_genes=extra_panel,
        containment=containment_for(families, best.locus) if families else (),
        quality=quality,
    )
