"""Unit tests for hit processing (SPEC.md §4, exact processing order)."""

import pytest

from gapit.blast import BlastRow
from gapit.hits import Hit, process_rows

TET_A_ID = "tinyamr~~~tetA~~~NC_000913.3:100-900~~~TETRACYCLINE"
TET_A_TITLE = TET_A_ID + " tetracycline efflux pump TetA"
SUL1_ID = "tinyamr~~~sul1~~~U12338.4:1-940~~~SULFONAMIDE"


def row(
    qseqid: str = "contig1",
    qstart: int = 1,
    qend: int = 100,
    qlen: int = 100,
    sseqid: str = TET_A_ID,
    sstart: int = 1,
    send: int = 100,
    slen: int = 100,
    sstrand: str = "plus",
    evalue: float = 1e-40,
    length: int = 100,
    pident: float = 100.0,
    gaps: int = 0,
    gapopen: int = 0,
    stitle: str = TET_A_TITLE,
) -> BlastRow:
    """A full-length 100%-identity tetA hit, with per-test overrides."""
    return BlastRow(
        qseqid=qseqid,
        qstart=qstart,
        qend=qend,
        qlen=qlen,
        sseqid=sseqid,
        sstart=sstart,
        send=send,
        slen=slen,
        sstrand=sstrand,
        evalue=evalue,
        length=length,
        pident=pident,
        gaps=gaps,
        gapopen=gapopen,
        stitle=stitle,
    )


def process(rows: list[BlastRow], mincov: float = 80.0, default_db: str = "tinyamr") -> list[Hit]:
    """process_rows with test-friendly defaults."""
    return process_rows(rows, mincov=mincov, default_db=default_db)


def test_full_length_hit_fields() -> None:
    """Given a canonical full-length hit, When processed, Then every field is
    populated from the row (coverage 100%, all-'=' map, '+' strand)."""
    (hit,) = process([row()])
    assert hit.sequence == "contig1"
    assert (hit.start, hit.end) == (1, 100)
    assert hit.strand == "+"
    assert hit.gene == "tetA"
    assert hit.database == "tinyamr"
    assert hit.accession == "NC_000913.3:100-900"
    assert hit.resistance == "TETRACYCLINE"
    assert hit.product == "tetracycline efflux pump TetA"
    assert (hit.s_start, hit.s_end, hit.s_len) == (1, 100, 100)
    assert hit.coverage_map == "==============="
    assert hit.gap_openings == 0
    assert hit.gaps == 0
    assert hit.identity_pct == 100.0
    assert hit.coverage_pct == 100.0


def test_minus_strand_swaps_subject_coords_only() -> None:
    """Given a minus-strand row with descending subject coords, When processed,
    Then s_start/s_end are swapped ascending, strand is '-', and the query
    coords are untouched."""
    hit = process([row(sstrand="minus", sstart=100, send=1, qstart=5, qend=104)])[0]
    assert (hit.s_start, hit.s_end) == (1, 100)
    assert (hit.start, hit.end) == (5, 104)
    assert hit.strand == "-"


def test_coverage_map_uses_post_swap_coords() -> None:
    """Given a full-length minus-strand hit, When processed, Then the map is
    computed from the ascending (post-swap) coords -> all '='."""
    hit = process([row(sstrand="minus", sstart=100, send=1)])[0]
    assert hit.coverage_map == "==============="


def test_dedup_first_row_wins_on_identical_query_span() -> None:
    """Given two genes aligned to the same (qseqid, qstart, qend), When
    processed, Then only the first BLAST row survives."""
    rows = [row(), row(sseqid=SUL1_ID)]
    hits = process(rows)
    assert len(hits) == 1
    assert hits[0].gene == "tetA"


def test_dedup_ignores_strand() -> None:
    """Given a plus-strand and a minus-strand row on the same query span, When
    processed, Then only the first survives with its own strand '+'."""
    rows = [row(), row(sstrand="minus", sstart=100, send=1)]
    hits = process(rows)
    assert len(hits) == 1
    assert hits[0].strand == "+"


def test_dedup_marks_key_even_when_first_row_is_filtered() -> None:
    """Given a below-threshold row followed by an above-threshold row on the
    same query span, When processed, Then zero hits — dedup (step 2) runs
    before the coverage filter (step 3) and the first row claims the key."""
    rows = [row(length=50), row()]
    assert process(rows) == []


def test_different_qstart_both_kept() -> None:
    """Given two hits on the same contig at different qstart, When processed,
    Then both survive dedup."""
    rows = [row(qstart=1), row(qstart=50, qend=100)]
    assert len(process(rows)) == 2


def test_coverage_exactly_threshold_is_kept() -> None:
    """Given 100*(80-0)/100 = 80.0%, When filtered at mincov=80, Then kept."""
    hits = process([row(length=80, qend=80)], mincov=80.0)
    assert len(hits) == 1
    assert hits[0].coverage_pct == 80.0


def test_coverage_just_below_threshold_is_dropped() -> None:
    """Given 100*(796-0)/1000 = 79.6%, When filtered at mincov=80, Then dropped."""
    assert process([row(length=796, slen=1000, qend=796, qlen=1000)]) == []


def test_coverage_79_996_displays_80_but_is_dropped() -> None:
    """Given 100*(20000-1)/25000 = 79.996% (displays as 80.00), When filtered
    at mincov=80, Then dropped — the comparison is on the unrounded float."""
    rows = [row(length=20000, gaps=1, slen=25000, qend=20000, qlen=20000)]
    assert process(rows) == []


def test_coverage_raw_float_is_stored_when_kept() -> None:
    """Given the 79.996% hit and mincov=79.99, When processed, Then it is kept
    with the raw unrounded coverage_pct."""
    rows = [row(length=20000, gaps=1, slen=25000, qend=20000, qlen=20000)]
    (hit,) = process(rows, mincov=79.99)
    assert hit.coverage_pct == pytest.approx(79.996)


def test_identity_is_never_filtered() -> None:
    """Given pident=10.0 with full coverage, When processed at default mincov,
    Then the hit is kept — identity is never post-filtered."""
    (hit,) = process([row(pident=10.0)])
    assert hit.identity_pct == 10.0


def test_product_strips_commas_and_tabs() -> None:
    """Given a stitle containing ',' and tab, When processed, Then both are
    stripped (no ~~~ present, so no token drop)."""
    hit = process([row(stitle="foo, bar\tbaz qux")])[0]
    assert hit.product == "foo barbaz qux"


def test_product_drops_leading_token_when_idsep_present() -> None:
    """Given a makeblastdb-style stitle prefixed with the ~~~ id, When
    processed, Then the leading whitespace-delimited token is dropped."""
    assert process([row()])[0].product == "tetracycline efflux pump TetA"


def test_product_empty_stitle_becomes_na() -> None:
    """Given an empty stitle, When processed, Then product is 'n/a'."""
    assert process([row(stitle="")])[0].product == "n/a"


def test_product_bare_idsep_token_without_description_stays() -> None:
    """Given a stitle that is only the ~~~ id (no trailing description), When
    processed, Then it is kept verbatim — the Perl leading-token regex requires
    trailing whitespace, so a bare id is not stripped (parity)."""
    assert process([row(stitle=TET_A_ID)])[0].product == TET_A_ID


def test_empty_sstrand_counts_as_plus() -> None:
    """Given sstrand='' (blastx rows), When processed, Then strand is '+'."""
    assert process([row(sstrand="")])[0].strand == "+"


def test_partial_header_falls_back_to_default_db() -> None:
    """Given a plain sseqid with no ~~~, When processed, Then gene is the whole
    id and database/accession/resistance come from the default_db fallback."""
    hit = process([row(sseqid="tetA(1)")], default_db="resfinder")[0]
    assert hit.gene == "tetA(1)"
    assert hit.database == "resfinder"
    assert hit.accession == ""
    assert hit.resistance == ""


def test_empty_rows_yield_no_hits() -> None:
    """Given no BLAST rows, When processed, Then no hits and no error."""
    assert process([]) == []


def test_stable_order_on_equal_sequence_start_keys() -> None:
    """Given two hits sharing (sequence, start), When Report-sorted, Then
    process order is preserved — the sort key must not extend to end/gene."""
    rows = [
        row(sseqid="db~~~geneA~~~acc~~~RES", qstart=5, qend=99),
        row(sseqid="db~~~geneB~~~acc~~~RES", qstart=5, qend=10),
    ]
    first, second = process(rows, mincov=0.0, default_db="db")

    def sort_key(hit: Hit) -> tuple[str, int]:
        return (hit.sequence, hit.start)

    assert sorted([first, second], key=sort_key) == [first, second]
    assert sorted([second, first], key=sort_key) == [second, first]
