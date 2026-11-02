"""Query/target alignment evidence for one cluster locus.

The coverage contract in :mod:`gapit.cluster_math` only needs target-space
runs.  The richer map (insertions, deletions, pieces, flanks, SVG) needs the
query side too, so this module keeps the full cs walk and projects it onto
the sample contig.  Coordinates are half-open: gene/feature coordinates stay
1-based inclusive at the model boundary, while ``target_position`` is the
0-based boundary between target bases.
"""

from collections.abc import Iterable
from typing import Literal

from pydantic import BaseModel

from gapit.cluster_math import parse_cs_steps
from gapit.cluster_profile import CoverageMode
from gapit.errors import GapitError
from gapit.gbfeatures import GeneFeature, LocusFeatures
from gapit.paf import PafRecord

Flank = Literal["between", "within_gene", "upstream", "downstream", "unknown"]


class AlignmentPiece(BaseModel, frozen=True):
    """One PAF row as a query/target block (forward contig coordinates)."""

    contig: str
    qstart: int
    qend: int
    tstart: int
    tend: int
    strand: Literal["+", "-"]
    matches: int
    alen: int
    identity_pct: float
    is_primary: bool


class Insertion(BaseModel, frozen=True):
    """A query-only segment (extra sample sequence) inside a locus."""

    contig: str
    qstart: int
    qend: int
    length: int
    target_position: int
    left_gene: str | None
    right_gene: str | None
    within_gene: str | None
    flank: Flank


class Deletion(BaseModel, frozen=True):
    """A target-only segment (reference sequence absent from the sample)."""

    contig: str
    qposition: int
    tstart: int
    tend: int
    length: int
    left_gene: str | None
    right_gene: str | None
    spanned_genes: tuple[str, ...]
    flank: Flank


class LocusAlignment(BaseModel, frozen=True):
    """All alignment evidence for one locus against one sample."""

    locus: str
    length: int
    pieces: tuple[AlignmentPiece, ...]
    insertions: tuple[Insertion, ...]
    deletions: tuple[Deletion, ...]


class QuerySpan(BaseModel, frozen=True):
    """A gene's detected span on a sample contig (forward coordinates)."""

    contig: str
    qstart: int
    qend: int


def _gene_interval(gene: GeneFeature) -> tuple[int, int]:
    """1-based inclusive feature -> 0-based half-open target interval."""
    return gene.start - 1, gene.end


def _insertion_flank(
    genes: tuple[GeneFeature, ...], position: int
) -> tuple[str | None, str | None, str | None, Flank]:
    """Nearest genes around a 0-based insertion boundary."""
    within = next(
        (g for g in genes if _gene_interval(g)[0] <= position < _gene_interval(g)[1]), None
    )
    left = max(
        (g for g in genes if _gene_interval(g)[1] <= position), key=lambda g: g.end, default=None
    )
    right = min(
        (g for g in genes if _gene_interval(g)[0] >= position), key=lambda g: g.start, default=None
    )
    if within is not None:
        return (
            left.gene_id if left else None,
            right.gene_id if right else None,
            within.gene_id,
            "within_gene",
        )
    if left is not None and right is not None:
        return left.gene_id, right.gene_id, None, "between"
    if right is not None:
        return None, right.gene_id, None, "upstream"
    if left is not None:
        return left.gene_id, None, None, "downstream"
    return None, None, None, "unknown"


def _deletion_flank(
    genes: tuple[GeneFeature, ...], start: int, end: int
) -> tuple[str | None, str | None, tuple[str, ...], Flank]:
    """Nearest genes and spanned genes for a 0-based deletion interval."""
    left = max(
        (g for g in genes if _gene_interval(g)[1] <= start), key=lambda g: g.end, default=None
    )
    right = min(
        (g for g in genes if _gene_interval(g)[0] >= end), key=lambda g: g.start, default=None
    )
    spanned = tuple(
        g.gene_id for g in genes if _gene_interval(g)[0] < end and _gene_interval(g)[1] > start
    )
    if spanned:
        flank: Flank = "within_gene"
    elif left is not None and right is not None:
        flank = "between"
    elif right is not None:
        flank = "upstream"
    elif left is not None:
        flank = "downstream"
    else:
        flank = "unknown"
    return (
        left.gene_id if left else None,
        right.gene_id if right else None,
        spanned,
        flank,
    )


def _forward_query_interval(row: PafRecord, qstart: int, qend: int) -> tuple[int, int]:
    """Alignment-orientation offsets -> forward contig half-open interval."""
    if row.strand == "+":
        return row.qstart + qstart, row.qstart + qend
    return row.qend - qend, row.qend - qstart


def _forward_query_position(row: PafRecord, qposition: int) -> int:
    """Alignment-orientation boundary -> forward contig boundary."""
    if row.strand == "+":
        return row.qstart + qposition
    return row.qend - qposition


def _piece(row: PafRecord) -> AlignmentPiece:
    identity = 100.0 * row.nmatch / row.alen if row.alen else 100.0
    return AlignmentPiece(
        contig=row.qname,
        qstart=row.qstart,
        qend=row.qend,
        tstart=row.tstart,
        tend=row.tend,
        strand=row.strand,
        matches=row.nmatch,
        alen=row.alen,
        identity_pct=identity,
        is_primary=row.is_primary,
    )


def _merge_spans(spans: list[QuerySpan]) -> tuple[QuerySpan, ...]:
    """Sort and merge overlapping/adjacent spans per contig."""
    merged: list[QuerySpan] = []
    for span in sorted(spans, key=lambda item: (item.contig, item.qstart, item.qend)):
        if merged and merged[-1].contig == span.contig and span.qstart <= merged[-1].qend:
            previous = merged[-1]
            merged[-1] = QuerySpan(
                contig=previous.contig,
                qstart=previous.qstart,
                qend=max(previous.qend, span.qend),
            )
        else:
            merged.append(span)
    return tuple(merged)


def gene_query_spans(
    locus: LocusFeatures,
    rows: Iterable[PafRecord],
    *,
    coverage_mode: CoverageMode = "primary",
) -> dict[str, tuple[QuerySpan, ...]]:
    """Query-side spans for every reference gene with aligned bases.

    Only match/substitution steps contribute: a reference gene is *present*
    when reference bases align, not when a query insertion happens to sit
    between its coordinates.  Spans are returned per gene, merged per contig.
    """
    collected: dict[str, list[QuerySpan]] = {}
    for row in rows:
        if row.tname != locus.id:
            continue
        if coverage_mode == "primary" and not row.is_primary:
            continue
        if row.cs is None:
            raise GapitError(
                f"minimap2 row without a cs tag (cannot map cluster {row.tname})",
                code="CS_MISSING",
                context={"locus": row.tname},
            )
        steps = parse_cs_steps(row.cs)
        for gene in locus.genes:
            gene_start, gene_end = _gene_interval(gene)
            for step in steps:
                if step.op not in ("match", "substitution"):
                    continue
                overlap_start = max(gene_start, step.tstart)
                overlap_end = min(gene_end, step.tend)
                if overlap_end <= overlap_start:
                    continue
                q_offset_start = step.qstart + (overlap_start - step.tstart)
                q_offset_end = step.qstart + (overlap_end - step.tstart)
                qstart, qend = _forward_query_interval(row, q_offset_start, q_offset_end)
                collected.setdefault(gene.gene_id, []).append(
                    QuerySpan(contig=row.qname, qstart=qstart, qend=qend)
                )
    return {gene_id: _merge_spans(spans) for gene_id, spans in collected.items()}


def build_locus_alignment(
    locus: LocusFeatures,
    rows: Iterable[PafRecord],
    *,
    coverage_mode: CoverageMode = "primary",
    min_indel_bp: int = 100,
) -> LocusAlignment:
    """Project every locus row's cs walk into pieces and indels.

    ``primary`` keeps the historical coverage row set; ``primary_secondary``
    includes secondary rows exactly like the coverage profile.  Indels below
    ``min_indel_bp`` are omitted from the map (they remain in the cs walk and
    in coverage math); the full walk is still available to callers that need
    every base.
    """
    pieces: list[AlignmentPiece] = []
    insertions: list[Insertion] = []
    deletions: list[Deletion] = []
    length = 0
    for row in rows:
        if row.tname != locus.id:
            continue
        if coverage_mode == "primary" and not row.is_primary:
            continue
        if not length:
            length = row.tlen
        if row.cs is None:
            raise GapitError(
                f"minimap2 row without a cs tag (cannot map cluster {row.tname})",
                code="CS_MISSING",
                context={"locus": row.tname},
            )
        pieces.append(_piece(row))
        for step in parse_cs_steps(row.cs):
            if step.op == "insertion":
                size = step.qend - step.qstart
                if size < min_indel_bp:
                    continue
                qstart, qend = _forward_query_interval(row, step.qstart, step.qend)
                left, right, within, flank = _insertion_flank(locus.genes, row.tstart + step.tstart)
                insertions.append(
                    Insertion(
                        contig=row.qname,
                        qstart=qstart,
                        qend=qend,
                        length=size,
                        target_position=row.tstart + step.tstart,
                        left_gene=left,
                        right_gene=right,
                        within_gene=within,
                        flank=flank,
                    )
                )
            elif step.op == "deletion":
                size = step.tend - step.tstart
                if size < min_indel_bp:
                    continue
                tstart = row.tstart + step.tstart
                tend = row.tstart + step.tend
                left, right, spanned, flank = _deletion_flank(locus.genes, tstart, tend)
                deletions.append(
                    Deletion(
                        contig=row.qname,
                        qposition=_forward_query_position(row, step.qstart),
                        tstart=tstart,
                        tend=tend,
                        length=size,
                        left_gene=left,
                        right_gene=right,
                        spanned_genes=spanned,
                        flank=flank,
                    )
                )
    pieces.sort(key=lambda piece: (piece.contig, piece.qstart, piece.tstart))
    insertions.sort(key=lambda item: (item.contig, item.target_position, item.qstart))
    deletions.sort(key=lambda item: (item.contig, item.tstart))
    return LocusAlignment(
        locus=locus.id,
        length=length,
        pieces=tuple(pieces),
        insertions=tuple(insertions),
        deletions=tuple(deletions),
    )
