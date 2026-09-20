"""Unit tests for PAF identity math and the gapit.reads/2 row filters.

Identity contract (SPEC.md §10, reads/2): per-alignment identity is
``100.0 * (alen - nm) / alen`` over PAF column 10 and the ``NM:i:`` tag.
A row without NM cannot be assessed and counts as 100.0; a degenerate
zero-length alignment block also counts as 100.0 (no division by zero).
"""

from gapit.paf import PafRecord, alignment_identity, filter_alignments, parse_paf_row


def paf(
    *,
    alen: int = 100,
    nm: int | None = 0,
    mapq: int = 60,
    tname: str = "db~~~geneA~~~ACC~~~RES",
) -> PafRecord:
    return PafRecord(
        qname="r1",
        qlen=alen,
        qstart=0,
        qend=alen,
        strand="+",
        tname=tname,
        tlen=100,
        tstart=0,
        tend=alen,
        nmatch=alen,
        alen=alen,
        mapq=mapq,
        nm=nm,
    )


# --- parsing -----------------------------------------------------------------


def test_parse_paf_row_extracts_nm_tag() -> None:
    """Given a PAF row with NM:i:7, When parsed, Then nm == 7."""
    line = (
        "r1\t100\t0\t50\t+\tdb~~~geneA~~~ACC~~~RES\t100\t10\t60\t50\t50\t60\t"
        "tp:A:P\tcm:i:18\tNM:i:7"
    )
    assert parse_paf_row(line).nm == 7


def test_parse_paf_row_without_nm_tag_is_none() -> None:
    """Given a PAF row with tags but no NM, When parsed, Then nm is None
    (minimap2 omits NM unless --cs is passed)."""
    line = (
        "r1\t100\t0\t50\t+\tdb~~~geneA~~~ACC~~~RES\t100\t10\t60\t50\t50\t60\t"
        "tp:A:P\tcm:i:18\tdv:f:0.0044"
    )
    assert parse_paf_row(line).nm is None


def test_parse_paf_row_bare_row_nm_is_none() -> None:
    """Given a 12-field row with no tags at all, When parsed, Then nm is None."""
    line = "r1\t100\t0\t50\t+\tdb~~~geneA~~~ACC~~~RES\t100\t10\t60\t50\t50\t60"
    record = parse_paf_row(line)
    assert record.nm is None
    assert record.is_primary is True


# --- identity arithmetic -------------------------------------------------------


def test_identity_is_mismatches_over_block_length() -> None:
    """Given alen=100 nm=5, When computed, Then identity is 95.0."""
    assert alignment_identity(paf(alen=100, nm=5)) == 95.0


def test_identity_perfect_alignment_is_100() -> None:
    assert alignment_identity(paf(alen=100, nm=0)) == 100.0


def test_identity_missing_nm_counts_as_100() -> None:
    """Given a row with no NM tag, When computed, Then identity is 100.0 —
    the row cannot be assessed and must not be filtered out."""
    assert alignment_identity(paf(alen=100, nm=None)) == 100.0


def test_identity_zero_alen_guard_counts_as_100() -> None:
    """Given alen=0 (degenerate block minimap2 never emits), When computed,
    Then identity is 100.0 instead of a ZeroDivisionError."""
    assert alignment_identity(paf(alen=0, nm=0)) == 100.0
    assert alignment_identity(paf(alen=0, nm=None)) == 100.0


# --- filtering -----------------------------------------------------------------


def test_filter_keeps_rows_at_or_above_identity_threshold() -> None:
    """Given rows at 95.0 and 94.99, When filtered at 95, Then only the
    exact-threshold row survives (>= comparison)."""
    kept = filter_alignments([paf(alen=100, nm=5), paf(alen=10000, nm=501)], min_identity=95.0)
    assert [row.nm for row in kept] == [5]


def test_filter_keeps_rows_at_or_above_mapq_threshold() -> None:
    """Given mapq 30 and 29, When filtered at 30, Then only the exact-threshold
    row survives."""
    kept = filter_alignments([paf(mapq=30), paf(mapq=29)], min_mapq=30)
    assert [row.mapq for row in kept] == [30]


def test_filter_combines_identity_and_mapq() -> None:
    """Given rows failing either predicate, When both filters are on, Then
    only rows passing both survive."""
    rows = [
        paf(alen=100, nm=5, mapq=30),
        paf(alen=100, nm=5, mapq=1),
        paf(alen=100, nm=50, mapq=60),
    ]
    kept = filter_alignments(rows, min_identity=95.0, min_mapq=30)
    assert len(kept) == 1
    assert kept[0].nm == 5 and kept[0].mapq == 30


def test_filter_off_keeps_every_row_unchanged() -> None:
    """Given both thresholds at their off defaults (0/0), When filtered,
    Then every row is kept verbatim — the gapit.reads/1 path must not drop a
    single row (byte-identity of the default output)."""
    rows = [paf(alen=100, nm=99, mapq=0), paf(alen=100, nm=None, mapq=0)]
    assert filter_alignments(rows) == rows


def test_filter_zero_identity_threshold_keeps_zero_identity_rows() -> None:
    """Given min_identity=0 (off) with min_mapq active, When filtered, Then a
    0%-identity row survives the identity side (>= 0)."""
    kept = filter_alignments([paf(alen=100, nm=100, mapq=5)], min_mapq=1)
    assert len(kept) == 1
