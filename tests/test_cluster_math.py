"""Pure-math tests for the cluster screening engine (TDD core of stage 2).

cs-tag parsing into target-space runs, per-gene coverage/identity over the
union of intersected alignment runs, gene verdicts, and the deterministic
locus ranking — all tested without minimap2 (the engine suites drive the
real binary).
"""

import pytest

from gapit.cluster_math import (
    PARTIAL_MIN_COV,
    gene_verdict,
    intersect_coverage,
    locus_sort_key,
    parse_cs,
)
from gapit.errors import GapitError


class TestParseCs:
    def test_match_run_only_when_single_exact_op(self) -> None:
        """Given a cs of one exact-match op, When parsed from tstart, Then it
        is one target run covering n bases with n matches."""
        assert parse_cs(":291", 0) == ((0, 291, 291),)

    def test_mixed_ops_walk_target_positions(self) -> None:
        """Given match/sub/ins/del ops, When parsed, Then runs walk target
        positions: sub consumes 1 target base with 0 matches, ins consumes 0,
        del consumes len bases with 0 matches."""
        runs = parse_cs(":10*ac+gt-cc:5", 100)
        assert runs == ((100, 110, 10), (110, 111, 0), (111, 113, 0), (113, 118, 5))

    def test_insertion_does_not_consume_target(self) -> None:
        """Given an insertion between two matches, When parsed, Then the
        target positions stay contiguous across it."""
        assert parse_cs(":4+acg:4", 10) == ((10, 14, 4), (14, 18, 4))

    def test_unknown_op_is_typed_error(self) -> None:
        """Given a cs with an op outside : * + -, When parsed, Then a
        CS_PARSE_FAILED error carries the offending op."""
        with pytest.raises(GapitError) as raised:
            parse_cs(":5=acgt", 0)
        assert raised.value.code == "CS_PARSE_FAILED"

    def test_malformed_length_is_typed_error(self) -> None:
        """Given a match op without an integer length, When parsed, Then a
        CS_PARSE_FAILED error names the op."""
        with pytest.raises(GapitError) as raised:
            parse_cs(":xy", 0)
        assert raised.value.code == "CS_PARSE_FAILED"

    def test_empty_cs_yields_no_runs(self) -> None:
        """Given an empty cs value, When parsed, Then no runs are produced."""
        assert parse_cs("", 0) == ()


class TestIntersectCoverage:
    def test_single_record_full_gene(self) -> None:
        """Given one record's runs covering the gene entirely with matches,
        When intersected, Then covered is the gene length and matched counts
        only the matching runs."""
        runs = ((0, 80, 80), (80, 90, 0), (90, 100, 10))
        assert intersect_coverage([runs], 0, 100) == (100, 90)

    def test_runs_clip_to_gene_interval(self) -> None:
        """Given runs spilling past the gene bounds, When intersected, Then
        only the gene slice counts (half-open clipping)."""
        runs = ((0, 150, 150),)
        assert intersect_coverage([runs], 50, 100) == (50, 50)

    def test_overlapping_records_dedup_per_base(self) -> None:
        """Given two records where the second mismatches a stretch the first
        matches, When intersected, Then the base is matched (any-record
        match wins) and identity is 100% over the union."""
        first = ((0, 60, 60),)
        second = ((50, 60, 0), (60, 100, 40))
        covered, matched = intersect_coverage([first, second], 0, 100)
        assert (covered, matched) == (100, 100)

    def test_cross_contig_fragments_union(self) -> None:
        """Given two records from different contigs covering complementary
        gene halves, When intersected, Then the union covers the whole gene
        (fragmentation resilience)."""
        left = ((0, 60, 60),)
        right = ((55, 100, 45),)
        assert intersect_coverage([left, right], 0, 100) == (100, 100)

    def test_no_runs_cover_zero(self) -> None:
        """Given records that never touch the gene, When intersected, Then
        (covered, matched) is (0, 0) — the identity-0 edge case."""
        runs = ((200, 300, 100),)
        assert intersect_coverage([runs], 0, 100) == (0, 0)


class TestVerdicts:
    def test_present_when_both_thresholds_met(self) -> None:
        """Given coverage and identity at/above their minimums, When called,
        Then the verdict is present (boundary-inclusive)."""
        assert gene_verdict(90.0, 90.0, min_gene_cov=90.0, min_gene_id=90.0) == "present"

    def test_partial_when_coverage_between_50_and_threshold(self) -> None:
        """Given 60% coverage of a gene needing 90, When called, Then the
        verdict is partial (the documented 50% floor)."""
        assert gene_verdict(60.0, 99.0, min_gene_cov=90.0, min_gene_id=90.0) == "partial"

    def test_partial_when_identity_below_threshold(self) -> None:
        """Given full coverage but identity under the minimum, When called,
        Then the verdict is partial, not present."""
        assert gene_verdict(100.0, 85.0, min_gene_cov=90.0, min_gene_id=90.0) == "partial"

    def test_absent_below_50_percent_coverage(self) -> None:
        """Given 49.9% coverage, When called, Then the verdict is absent."""
        assert gene_verdict(49.9, 100.0, min_gene_cov=90.0, min_gene_id=90.0) == "absent"
        assert PARTIAL_MIN_COV == 50.0


class TestLocusRanking:
    def test_coverage_desc_dominates(self) -> None:
        """Given two loci, When ranked, Then higher coverage sorts first
        regardless of identity."""
        key = locus_sort_key
        assert key("locusB", 99.0, 99.0, 900) < key("locusA", 98.0, 100.0, 1200)

    def test_identity_breaks_coverage_ties(self) -> None:
        """Given equal coverage, When ranked, Then higher identity sorts
        first."""
        key = locus_sort_key
        assert key("locusB", 99.0, 98.0, 900) < key("locusA", 99.0, 97.0, 1200)

    def test_covered_bp_breaks_identity_ties(self) -> None:
        """Given equal coverage and identity, When ranked, Then more covered
        bases sort first (short-reference bias guard)."""
        key = locus_sort_key
        assert key("locusB", 99.0, 99.0, 1200) < key("locusA", 99.0, 99.0, 900)

    def test_locus_id_breaks_remaining_ties_deterministically(self) -> None:
        """Given identical stats, When ranked, Then the lexically smaller
        locus id sorts first."""
        key = locus_sort_key
        assert key("locusA", 99.0, 99.0, 900) < key("locusB", 99.0, 99.0, 900)
