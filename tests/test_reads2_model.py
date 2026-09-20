"""Unit + model tests for gapit.reads/2: aggregation weighting, document
shape, and the byte-identity of gapit.reads/1 rendering when the new
thresholds are off.
"""

import json
from datetime import UTC, datetime

from pydantic import TypeAdapter

from gapit.formats.json import (
    Reads2Document,
    render_reads2_json,
    render_reads_json,
)
from gapit.paf import PafRecord
from gapit.reads import GeneCoverage, ReadsParams, ReadsReport, aggregate_coverage

PINNED_NOW = datetime(2026, 9, 20, 12, 0, 0, tzinfo=UTC)
GENE = "db~~~geneA~~~ACC~~~RES"


def paf_row(
    alen: int, nm: int | None, *, tstart: int = 0, tend: int = 50, qname: str = "r1"
) -> PafRecord:
    return PafRecord(
        qname=qname,
        qlen=alen,
        qstart=0,
        qend=alen,
        strand="+",
        tname=GENE,
        tlen=100,
        tstart=tstart,
        tend=tend,
        nmatch=alen,
        alen=alen,
        mapq=60,
        nm=nm,
    )


# --- aggregation weighting ---------------------------------------------------


def test_mean_identity_is_alen_weighted() -> None:
    """Given primary rows at 100% (alen 100) and 90% (alen 300), When
    aggregated, Then mean_identity_pct is the alen-weighted mean 92.5."""
    (entry,) = aggregate_coverage(
        [paf_row(100, 0), paf_row(300, 30, tstart=50, tend=100, qname="r2")],
        default_db="db",
        min_breadth=90.0,
    )
    assert entry.mean_identity_pct == 92.5


def test_mean_identity_counts_missing_nm_rows_as_100() -> None:
    """Given one row without NM (100.0) and one at 50%, When aggregated,
    Then the weighted mean is 75.0."""
    (entry,) = aggregate_coverage(
        [paf_row(100, None), paf_row(100, 50, tstart=50, tend=100, qname="r2")],
        default_db="db",
        min_breadth=90.0,
    )
    assert entry.mean_identity_pct == 75.0


def test_mean_identity_excludes_secondary_rows() -> None:
    """Given a secondary low-identity row, When aggregated, Then only the
    primary row feeds the mean."""
    rows = [
        paf_row(100, 0),
        PafRecord(
            qname="r2",
            qlen=100,
            qstart=0,
            qend=100,
            strand="+",
            tname=GENE,
            tlen=100,
            tstart=0,
            tend=100,
            nmatch=100,
            alen=100,
            mapq=0,
            nm=90,
            is_primary=False,
        ),
    ]
    (entry,) = aggregate_coverage(rows, default_db="db", min_breadth=90.0)
    assert entry.mean_identity_pct == 100.0


# --- document shape -----------------------------------------------------------


def sample_report() -> tuple[GeneCoverage, ReadsParams]:
    gene = GeneCoverage(
        database="db",
        gene="geneA",
        accession="ACC",
        function="RES",
        product="demo product",
        tlen=100,
        breadth_pct=97.5,
        mean_depth=1.2,
        reads_mapped=3,
        present=True,
        mean_identity_pct=98.756,
    )
    params = ReadsParams(
        db="db", read_type="sr", min_breadth=90.0, threads=1, min_identity=95.0, min_mapq=0
    )
    return gene, params


def test_reads2_json_document_shape() -> None:
    """Given a filtered run rendered as reads/2, When parsed, Then the schema
    literal is gapit.reads/2, params carry both thresholds, and gene entries
    carry mean_identity_pct rounded to 2."""
    gene, params = sample_report()
    report = ReadsReport(reads=("a.fq",), genes=(gene,))
    output = render_reads2_json([report], params, now=PINNED_NOW)
    document = TypeAdapter(Reads2Document).validate_json(output)
    assert document.schema_name == "gapit.reads/2"
    assert document.params.min_identity == 95.0
    assert document.params.min_mapq == 0
    assert document.params.db == "db"
    (entry,) = document.files[0].genes
    assert entry.mean_identity_pct == 98.76
    assert entry.gene == "geneA"
    assert entry.present is True
    assert document.files[0].reads == ["a.fq"]


def test_reads1_rendering_omits_reads2_fields() -> None:
    """Given params carrying off-thresholds, When rendered as reads/1, Then
    the JSON has no min_identity/min_mapq/mean_identity_pct keys at all."""
    gene, _ = sample_report()
    params = ReadsParams(db="db")
    report = ReadsReport(reads=("a.fq",), genes=(gene,))
    output = render_reads_json([report], params, now=PINNED_NOW)
    assert '"schema": "gapit.reads/1"' in output
    assert "min_identity" not in output
    assert "min_mapq" not in output
    assert "mean_identity_pct" not in output


def test_reads1_rendering_ignores_nonzero_thresholds() -> None:
    """Given reads/1 rendering of params that (hypothetically) carry nonzero
    thresholds, When rendered, Then the /1 document shows only the frozen /1
    params (rendering selects fields; selection of /1 itself never happens
    with thresholds on)."""
    gene, params = sample_report()
    report = ReadsReport(reads=("a.fq",), genes=(gene,))
    output = render_reads_json([report], params, now=PINNED_NOW)
    document = json.loads(output)
    assert set(document["params"]) == {"db", "read_type", "min_breadth", "threads"}
