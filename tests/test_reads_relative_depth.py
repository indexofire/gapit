"""Relative-depth presence gating and the DEPTH_RATIO table column."""

from gapit.formats.reads_tsv import reads_tsv_chunk, reads_tsv_preamble
from gapit.paf import PafRecord
from gapit.reads import GeneCoverage, ReadsReport, aggregate_coverage

LPXA = "db~~~lpxA~~~ACC~~~CONTROL"
GENE_A = "db~~~geneA~~~ACC~~~RES"
GENE_B = "db~~~geneB~~~ACC~~~RES"


def row(tname: str, *, start: int, end: int, tlen: int = 100, qname: str = "r") -> PafRecord:
    return PafRecord(
        qname=qname,
        qlen=end - start,
        qstart=0,
        qend=end - start,
        strand="+",
        tname=tname,
        tlen=tlen,
        tstart=start,
        tend=end,
        nmatch=end - start,
        alen=end - start,
        mapq=60,
    )


def aggregate(
    rows: list[PafRecord],
    *,
    min_depth: float = 0.0,
    min_relative_depth: float = 0.0,
    relative_to: str | None = None,
) -> list[GeneCoverage]:
    return aggregate_coverage(
        rows,
        default_db="db",
        min_breadth=90.0,
        min_depth=min_depth,
        min_relative_depth=min_relative_depth,
        relative_to=relative_to,
    )


def test_relative_depth_ratio_gates_presence() -> None:
    """Given control depth 0.5, geneA 0.5 (100%) and geneB 0.1 (20%), When the
    threshold is 20%, Then geneA passes and geneB is exactly at the boundary
    (>= semantics) while a 10% gene fails."""
    rows = [row(LPXA, start=0, end=90, qname=f"c{i}") for i in range(10)]
    rows += [row(GENE_A, start=0, end=90, qname=f"a{i}") for i in range(10)]
    rows += [row(GENE_B, start=0, end=90, qname=f"b{i}") for i in range(2)]
    entries = aggregate(rows, min_relative_depth=0.20, relative_to="lpxA")
    by_gene = {entry.gene: entry for entry in entries}
    assert by_gene["lpxA"].relative_depth_pct == 100.0
    assert by_gene["geneA"].relative_depth_pct == 100.0
    assert by_gene["geneA"].present is True
    assert by_gene["geneB"].relative_depth_pct == 20.0
    assert by_gene["geneB"].present is True


def test_relative_depth_below_threshold_is_absent() -> None:
    rows = [row(LPXA, start=0, end=90, qname=f"c{i}") for i in range(10)]
    rows += [row(GENE_A, start=0, end=90, qname="a")]
    (gene,) = [
        e for e in aggregate(rows, min_relative_depth=0.20, relative_to="lpxA") if e.gene == "geneA"
    ]
    assert gene.relative_depth_pct == 10.0
    assert gene.present is False


def test_missing_control_makes_everything_absent() -> None:
    """Given no lpxA control, When relative gating is on, Then no call passes:
    the control's absence must not silently pass every gene."""
    (gene,) = aggregate([row(GENE_A, start=0, end=90)], min_relative_depth=0.20, relative_to="lpxA")
    assert gene.relative_depth_pct is None
    assert gene.present is False


def test_absolute_min_depth_is_an_independent_floor() -> None:
    rows = [row(LPXA, start=0, end=90), row(GENE_A, start=0, end=90, qname="a")]
    (gene,) = [e for e in aggregate(rows, min_depth=1.0, relative_to="lpxA") if e.gene == "geneA"]
    assert gene.mean_depth == 0.9
    assert gene.present is False


def test_tsv_appends_depth_ratio_in_relative_mode() -> None:
    report = ReadsReport(
        reads=("sample",),
        genes=(
            GeneCoverage(
                database="db",
                gene="geneA",
                accession="ACC",
                function="RES",
                product="product",
                tlen=100,
                breadth_pct=100.0,
                mean_depth=12.0,
                reads_mapped=6,
                present=True,
                mean_identity_pct=99.0,
                relative_depth_pct=40.0,
            ),
        ),
    )
    header = reads_tsv_preamble(reads2=True, relative=True)
    assert header.rstrip("\n").endswith("IDENTITY%\tDEPTH_RATIO")
    line = reads_tsv_chunk(report, reads2=True, relative=True)
    assert line.rstrip("\n").split("\t")[-2:] == ["99.00", "40.00"]


def test_reads2_json_carries_relative_depth_when_active() -> None:
    from datetime import UTC, datetime

    from gapit.formats.reads_json import render_reads2_json
    from gapit.reads import ReadsParams

    report = ReadsReport(
        reads=("sample",),
        genes=(
            GeneCoverage(
                database="db",
                gene="geneA",
                accession="ACC",
                function="RES",
                product="product",
                tlen=100,
                breadth_pct=100.0,
                mean_depth=12.0,
                reads_mapped=6,
                present=True,
                mean_identity_pct=99.0,
                relative_depth_pct=40.0,
            ),
        ),
    )
    params = ReadsParams(
        db="db",
        min_identity=80.0,
        min_mapq=20,
        min_depth=5.0,
        min_relative_depth=0.20,
        relative_to="lpxA",
    )
    document = render_reads2_json([report], params, now=datetime(2026, 10, 29, tzinfo=UTC))
    assert '"relative_depth_pct": 40.0' in document
    assert '"min_relative_depth": 0.2' in document
    assert '"relative_to": "lpxA"' in document
