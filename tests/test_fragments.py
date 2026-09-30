"""Unit tests for cross-contig gene-fragment merging (gapit extension).

The merge math is pure: interval unions over subject coordinates, the
minimap.py 15-char binning applied to a union of spans, aligned-length
weighted identity, and the anchor rule for query coordinates.
"""

from typing import Literal

import pytest

from gapit.fragments import (
    merge_fragment_hits,
    union_coverage_map,
    union_coverage_pct,
    union_spans,
)
from gapit.hits import Hit
from gapit.minimap import minimap


def frag(
    contig: str,
    qstart: int,
    qend: int,
    s_start: int,
    s_end: int,
    *,
    strand: Literal["+", "-"] = "+",
    pident: float = 100.0,
    slen: int = 600,
    gaps: int = 0,
    gapopen: int = 0,
    aligned: int | None = None,
    gene: str = "demoGene",
) -> Hit:
    """A decoded partial hit with per-test overrides; ``aligned`` defaults to
    the gap-free subject span length (blast length == span when ungapped)."""
    aligned_len = s_end - s_start + 1 if aligned is None else aligned
    return Hit(
        sequence=contig,
        start=qstart,
        end=qend,
        strand=strand,
        gene=gene,
        database="fragdb",
        accession="FRAG-0001",
        product="demo gene",
        function="DEMOCLASS",
        s_start=s_start,
        s_end=s_end,
        s_len=slen,
        coverage_map="===============.",
        gap_openings=gapopen,
        gaps=gaps,
        identity_pct=pident,
        coverage_pct=100.0 * aligned_len / slen,
        aligned_len=aligned_len,
    )


# --- interval union -------------------------------------------------------


def test_union_spans_merges_overlaps_not_abutting_pairs() -> None:
    """Given overlapping, contained, abutting, and unsorted spans, When
    unioned, Then overlaps collapse, abutting pairs stay separate (length is
    unaffected either way), and output is start-sorted."""
    assert union_spans([(10, 20), (1, 5)]) == [(1, 5), (10, 20)]
    assert union_spans([(1, 10), (10, 20)]) == [(1, 20)]
    assert union_spans([(1, 10), (11, 20)]) == [(1, 10), (11, 20)]
    assert union_spans([(1, 10), (2, 5)]) == [(1, 10)]


def test_union_coverage_pct_counts_each_subject_position_once() -> None:
    """Given disjoint (240 + 360) and overlapping (300 + 401) 600nt spans,
    When coverage is computed, Then it is 100*union_len/slen, unrounded."""
    assert union_coverage_pct([(1, 240), (241, 600)], 600) == 100.0
    assert union_coverage_pct([(1, 300), (260, 600)], 600) == 100.0
    assert union_coverage_pct([(1, 240), (241, 479)], 600) == pytest.approx(79.8333333)


def test_union_coverage_pct_boundary_is_inclusive() -> None:
    """Given a union of exactly 80% of the subject, When compared to
    mincov=80, Then the merge keeps it (>= on the unrounded float)."""
    assert union_coverage_pct([(1, 240), (241, 480)], 600) == 80.0


# --- union coverage map ---------------------------------------------------


def test_union_map_single_span_matches_minimap_exactly() -> None:
    """Given single spans, When rendered as a union map, Then the result is
    byte-identical to minimap.py for the same span (incl. int() quirks and
    the broken '/' variant)."""
    for span, broken in (
        ((1, 600), 0),
        ((10, 300), 0),
        ((2, 47), 0),
        ((10, 300), 2),
        ((1, 600), 1),
    ):
        assert union_coverage_map([span], 600, broken=broken > 0) == minimap(
            span[0], span[1], 600, broken
        )


def test_union_map_ors_disjoint_spans_into_15_chars() -> None:
    """Given disjoint spans [(1,180),(361,600)] on a 600nt subject, When
    rendered, Then covered bins OR together into exactly 15 chars with the
    uncovered middle as dots (scale 40: bins 0-4 and 9-15)."""
    assert union_coverage_map([(1, 180), (361, 600)], 600, broken=False) == "=====....======"


def test_union_map_overlapping_spans_mark_bin_once() -> None:
    """Given overlapping spans in the same bins, When rendered, Then the map
    shows '=' once per bin (no doubling) and stays 15 chars."""
    assert union_coverage_map([(1, 120), (100, 200)], 600, broken=False) == "======........."


def test_union_map_broken_variant_inserts_slash_at_half_width() -> None:
    """Given a union with any gap openings, When rendered, Then the map uses
    14 bins with '/' after bin 7 — the minimap broken-map rule applied to the
    union."""
    coverage_map = union_coverage_map([(1, 180), (361, 600)], 600, broken=True)
    assert len(coverage_map) == 15
    assert coverage_map.count("/") == 1
    assert coverage_map.index("/") == 8


# --- merge_fragment_hits --------------------------------------------------


def test_disjoint_fragments_merge_into_one_hit() -> None:
    """Given 40% and 60% partials of demoGene on two contigs (union 100%),
    When merged at mincov=80, Then ONE merged hit with union coverage, the
    comma-joined contig names, and both fragments attached."""
    partials = [frag("contigB", 1, 360, 241, 600), frag("contigA", 1, 240, 1, 240)]
    (merged,) = merge_fragment_hits([], partials, mincov=80.0)
    assert merged.merged is True
    assert merged.gene == "demoGene"
    assert merged.coverage_pct == 100.0
    assert merged.sequence == "contigA,contigB"
    assert [(f.contig, f.start, f.end) for f in merged.fragments] == [
        ("contigA", 1, 240),
        ("contigB", 1, 360),
    ]
    assert merged.s_start == 1
    assert merged.s_end == 600
    assert merged.coverage_map == "==============="


def test_union_below_mincov_reports_nothing() -> None:
    """Given 30% and 40% partials (union 70%), When merged at mincov=80,
    Then the gene stays unreported."""
    partials = [frag("contigA", 1, 180, 1, 180), frag("contigB", 1, 240, 361, 600)]
    assert merge_fragment_hits([], partials, mincov=80.0) == []


def test_union_exactly_at_mincov_is_kept() -> None:
    """Given two 40% partials whose union is exactly 80.0%, When merged at
    mincov=80, Then the merged hit is emitted (>= boundary)."""
    partials = [frag("contigA", 1, 240, 1, 240), frag("contigB", 1, 240, 241, 480)]
    (merged,) = merge_fragment_hits([], partials, mincov=80.0)
    assert merged.coverage_pct == 80.0


def test_weighted_identity_uses_aligned_lengths() -> None:
    """Given fragments at 90%/95% identity with aligned lengths 240/360,
    When merged, Then identity_pct is the aligned-length-weighted mean
    (90*240+95*360)/600 = 93.0."""
    partials = [
        frag("contigA", 1, 240, 1, 240, pident=90.0, aligned=240),
        frag("contigB", 1, 360, 241, 600, pident=95.0, aligned=360),
    ]
    (merged,) = merge_fragment_hits([], partials, mincov=80.0)
    assert merged.identity_pct == pytest.approx(93.0)


def test_anchor_is_the_largest_aligned_fragment() -> None:
    """Given a 240nt and a 360nt fragment (minus strand on the larger), When
    merged, Then START/END/STRAND come from the largest- aligned fragment —
    the documented anchor rule."""
    partials = [
        frag("contigA", 1, 240, 1, 240, strand="+"),
        frag("contigB", 5, 364, 241, 600, strand="-"),
    ]
    (merged,) = merge_fragment_hits([], partials, mincov=80.0)
    assert (merged.start, merged.end, merged.strand) == (5, 364, "-")


def test_anchor_tie_breaks_by_contig_then_start() -> None:
    """Given two equal-sized fragments, When merged, Then the anchor is the
    first in stable (contig, start) order — determinism."""
    partials = [
        frag("contigB", 1, 300, 301, 600, strand="-"),
        frag("contigA", 1, 300, 1, 300, strand="+"),
    ]
    (merged,) = merge_fragment_hits([], partials, mincov=80.0)
    assert (merged.start, merged.end, merged.strand) == (1, 300, "+")


def test_minus_strand_fragments_union_on_swapped_coords() -> None:
    """Given two minus-strand partials whose subject intervals were already
    swapped ascending (SPEC §4 step 1), When merged, Then the union is over
    the ascending intervals and the anchor strand is '-'."""
    partials = [
        frag("contigA", 1, 300, 1, 300, strand="-"),
        frag("contigB", 1, 300, 301, 600, strand="-"),
    ]
    (merged,) = merge_fragment_hits([], partials, mincov=80.0)
    assert merged.coverage_pct == 100.0
    assert merged.strand == "-"


def test_gene_with_full_hit_passes_through_untouched() -> None:
    """Given a gene that already has one >=mincov hit plus a stray partial,
    When merged, Then the normal hit is returned unchanged and the partial is
    NOT folded in (merging only rescues unreported genes)."""
    full = frag("contigA", 1, 600, 1, 600)
    partial = frag("contigB", 1, 100, 500, 599)
    result = merge_fragment_hits([full], [partial], mincov=80.0)
    assert result == [full]


def test_two_genes_merge_independently() -> None:
    """Given partials of two genes, When merged, Then each gene yields its own
    merged hit (grouping is per gene)."""
    partials = [
        frag("c1", 1, 240, 1, 240, gene="geneA"),
        frag("c2", 1, 360, 241, 600, gene="geneA"),
        frag("c3", 1, 200, 1, 200, gene="geneB"),
        frag("c4", 1, 400, 201, 600, gene="geneB"),
    ]
    merged = merge_fragment_hits([], partials, mincov=80.0)
    assert [hit.gene for hit in merged] == ["geneA", "geneB"]
    assert all(hit.merged for hit in merged)


def test_full_hits_keep_their_order_and_merge_results_are_appended() -> None:
    """Given two full hits and one mergeable gene, When merged, Then the full
    hits stay in input order and the merged hit follows (final ordering is the
    caller's (sequence, start) sort)."""
    full_a = frag("cA", 1, 600, 1, 600, gene="geneA")
    full_b = frag("cB", 1, 600, 1, 600, gene="geneB")
    partials = [
        frag("c1", 1, 240, 1, 240, gene="geneC"),
        frag("c2", 1, 360, 241, 600, gene="geneC"),
    ]
    result = merge_fragment_hits([full_a, full_b], partials, mincov=80.0)
    assert result[:2] == [full_a, full_b]
    assert result[2].gene == "geneC"


def test_gaps_and_gap_openings_are_summed() -> None:
    """Given fragments with gaps 2+3 and gap openings 1+2, When merged, Then
    GAPS totals are 3/5 and the union map takes the broken (14-box + '/') form."""
    partials = [
        frag("contigA", 1, 240, 1, 240, gaps=2, gapopen=1),
        frag("contigB", 1, 360, 241, 600, gaps=3, gapopen=2),
    ]
    (merged,) = merge_fragment_hits([], partials, mincov=80.0)
    assert (merged.gap_openings, merged.gaps) == (3, 5)
    assert "/" in merged.coverage_map
    assert len(merged.coverage_map) == 15


def test_fragments_are_sorted_by_contig_then_start() -> None:
    """Given partials in reverse order, When merged, Then the fragments array
    is sorted by (contig, start) regardless of input order."""
    partials = [
        frag("contigZ", 1, 360, 241, 600),
        frag("contigA", 50, 289, 1, 240),
        frag("contigA", 1, 49, 1, 240),
    ]
    (merged,) = merge_fragment_hits([], partials, mincov=80.0)
    assert [(f.contig, f.start) for f in merged.fragments] == [
        ("contigA", 1),
        ("contigA", 50),
        ("contigZ", 1),
    ]


def test_single_partial_cannot_merge() -> None:
    """Given a lone sub-threshold partial, When merged, Then nothing is
    reported — a group needs >= 2 fragments by contract."""
    assert merge_fragment_hits([], [frag("c1", 1, 240, 1, 240)], mincov=80.0) == []


def test_no_partials_returns_hits_unchanged() -> None:
    """Given only full hits, When merged, Then the list is returned as-is."""
    full = frag("cA", 1, 600, 1, 600)
    assert merge_fragment_hits([full], [], mincov=80.0) == [full]
