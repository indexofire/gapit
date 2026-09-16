"""Unit tests for per-gene coverage aggregation (SPEC.md §10).

All breadth/depth expectations are computed BY HAND from interval unions.
"""

from gapit.reads import GeneCoverage, PafRecord, aggregate_coverage

GENE_A = "db~~~geneA~~~ACC~~~RES"
GENE_B = "db~~~geneB~~~ACC~~~RES"


def paf(
    qname: str = "r1",
    tname: str = GENE_A,
    tlen: int = 100,
    tstart: int = 0,
    tend: int = 50,
    is_primary: bool = True,
) -> PafRecord:
    """A primary alignment of `qname` on `tname` spanning [tstart, tend)."""
    return PafRecord(
        qname=qname,
        qlen=tend - tstart,
        qstart=0,
        qend=tend - tstart,
        strand="+",
        tname=tname,
        tlen=tlen,
        tstart=tstart,
        tend=tend,
        nmatch=tend - tstart,
        alen=tend - tstart,
        mapq=60,
        is_primary=is_primary,
    )


def aggregate(rows: list[PafRecord], min_breadth: float = 90.0) -> list[GeneCoverage]:
    return aggregate_coverage(rows, default_db="db", min_breadth=min_breadth)


def test_overlapping_intervals_union_breadth() -> None:
    """Given [0,50) + [25,75) on tlen 100, When aggregated, Then breadth is
    75% (union) and mean depth is (50+50)/100 = 1.0."""
    (entry,) = aggregate([paf(tstart=0, tend=50), paf(qname="r2", tstart=25, tend=75)])
    assert entry.gene == "geneA"
    assert entry.breadth_pct == 75.0
    assert entry.mean_depth == 1.0
    assert entry.tlen == 100


def test_nested_interval_adds_depth_not_breadth() -> None:
    """Given [0,50) + [10,20) on tlen 100, When aggregated, Then breadth stays
    50% and depth becomes (50+10)/100 = 0.6."""
    (entry,) = aggregate([paf(tstart=0, tend=50), paf(qname="r2", tstart=10, tend=20)])
    assert entry.breadth_pct == 50.0
    assert entry.mean_depth == 0.6


def test_adjacent_intervals_cover_whole_target() -> None:
    """Given [0,50) + [50,100), When aggregated, Then breadth is 100%
    (adjacent half-open intervals tile the target)."""
    (entry,) = aggregate([paf(tstart=0, tend=50), paf(qname="r2", tstart=50, tend=100)])
    assert entry.breadth_pct == 100.0
    assert entry.present is True


def test_both_mates_counted_as_distinct_reads() -> None:
    """Given two qnames on the same span, When aggregated, Then reads_mapped
    is 2; a repeated qname counts once."""
    (entry,) = aggregate(
        [paf(qname="mate1"), paf(qname="mate2"), paf(qname="mate1", tstart=60, tend=80)]
    )
    assert entry.reads_mapped == 2


def test_multi_target_rows_aggregate_independently() -> None:
    """Given rows on geneA (75%) and geneB (50%), When aggregated, Then two
    entries with independent stats, sorted by breadth descending."""
    entries = aggregate(
        [
            paf(tstart=0, tend=50),
            paf(qname="r2", tstart=25, tend=75),
            paf(qname="r3", tname=GENE_B, tstart=0, tend=50),
        ]
    )
    assert [(e.gene, e.breadth_pct) for e in entries] == [("geneA", 75.0), ("geneB", 50.0)]


def test_sort_ties_break_on_gene_name() -> None:
    """Given equal breadth on geneB and geneA, When aggregated, Then the
    alphabetical gene order wins the tie."""
    entries = aggregate([paf(tname=GENE_B), paf(qname="r2", tname=GENE_A)])
    assert [e.gene for e in entries] == ["geneA", "geneB"]


def test_breadth_exactly_at_threshold_is_present() -> None:
    """Given 90/100 covered, When min_breadth=90, Then present (>= comparison)."""
    (entry,) = aggregate([paf(tstart=0, tend=90)])
    assert entry.breadth_pct == 90.0
    assert entry.present is True


def test_breadth_just_below_threshold_is_absent() -> None:
    """Given 8999/10000 covered (89.99%), When min_breadth=90, Then absent."""
    (entry,) = aggregate([paf(tlen=10000, tstart=0, tend=8999)])
    assert entry.breadth_pct == 89.99
    assert entry.present is False


def test_secondary_only_gene_is_omitted() -> None:
    """Given only non-primary rows for geneB, When aggregated, Then geneB is
    omitted entirely (zero-read genes are not reported)."""
    entries = aggregate([paf(), paf(qname="r2", tname=GENE_B, is_primary=False)])
    assert [e.gene for e in entries] == ["geneA"]


def test_zero_rows_yield_no_genes() -> None:
    assert aggregate([]) == []


def test_plain_tname_falls_back_to_default_db() -> None:
    """Given a tname without ~~~, When aggregated, Then gene is the whole id
    and database is the default."""
    (entry,) = aggregate([paf(tname="tetA(1)")])
    assert entry.gene == "tetA(1)"
    assert entry.database == "db"
    assert entry.accession == ""
