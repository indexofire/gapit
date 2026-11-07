"""Rich LocusMap assembly and JSON surface tests."""

from __future__ import annotations

from collections.abc import Iterator
from datetime import UTC, datetime
from pathlib import Path

import pytest

from gspec.cluster import (
    ClusterParams,
    ClusterReport,
    GeneCall,
    LocusCall,
    screen_cluster_file,
)
from gspec.cluster_alignment import QuerySpan, build_locus_alignment
from gspec.db import Database
from gspec.formats.locus_map import render_locus_map_json
from gspec.gbfeatures import FeaturesDocument, GeneFeature, LocusFeatures
from gspec.gene_families import GeneFamily, GeneFamilyDocument, LocusFamilySet
from gspec.locus_map import build_locus_map
from gspec.paf import PafRecord

LOCUS = LocusFeatures(
    id="OL1",
    label="O1",
    type="O1",
    genes=(
        GeneFeature(gene_id="g1", start=1, end=100, strand="+", product="gene 1"),
        GeneFeature(gene_id="g2", start=201, end=300, strand="+", product="gene 2"),
        GeneFeature(gene_id="g3", start=401, end=500, strand="+", product="gene 3"),
    ),
)
FEATURES = FeaturesDocument(loci=(LOCUS,))


def _call(*, missing: tuple[str, ...] = ()) -> LocusCall:
    genes = (
        GeneCall(
            gene_id="g1",
            start=1,
            end=100,
            strand="+",
            coverage_pct=100.0,
            identity_pct=100.0,
            verdict="absent" if "g1" in missing else "present",
        ),
        GeneCall(
            gene_id="g2",
            start=201,
            end=300,
            strand="+",
            coverage_pct=100.0,
            identity_pct=100.0,
            verdict="absent" if "g2" in missing else "present",
        ),
        GeneCall(
            gene_id="g3",
            start=401,
            end=500,
            strand="+",
            coverage_pct=100.0,
            identity_pct=100.0,
            verdict="absent" if "g3" in missing else "present",
        ),
    )
    return LocusCall(
        locus="OL1",
        label="O1",
        type="O1",
        coverage_pct=99.0,
        identity_pct=99.0,
        rank=1,
        genes=genes,
        missing=missing,
    )


def _row(cs: str) -> PafRecord:
    return PafRecord(
        qname="ctg1",
        qlen=1200,
        qstart=100,
        qend=850,
        strand="+",
        tname="OL1",
        tlen=500,
        tstart=0,
        tend=500,
        nmatch=750,
        alen=750,
        mapq=60,
        cs=cs,
    )


def test_locus_map_carries_genes_insertions_and_quality() -> None:
    """Given a high-quality locus with one insertion, When mapped, Then the
    map carries gene evidence, the insertion, and a passing quality call."""
    alignment = build_locus_alignment(LOCUS, [_row(f":100+{'A' * 150}:200+{'C' * 100}:200")])
    mapped = build_locus_map(_call(), alignment, FEATURES, (_call(),))
    assert mapped.length == 500
    assert [gene.product for gene in mapped.genes] == ["gene 1", "gene 2", "gene 3"]
    assert len(mapped.insertions) == 2
    assert mapped.quality.status == "pass"
    assert mapped.quality.extra_bp == 250
    assert mapped.quality.missing_genes == 0


def test_locus_map_quality_fails_on_missing_genes() -> None:
    """Given more missing genes than the prototype allowance, When mapped,
    Then quality fails with a reason but the phenotype surface is untouched."""
    alignment = build_locus_alignment(LOCUS, [_row(":500")])
    mapped = build_locus_map(
        _call(missing=("g1", "g2")),
        alignment,
        FEATURES,
        (_call(missing=("g1", "g2")),),
        max_missing_genes=1,
    )
    assert mapped.quality.status == "fail"
    assert "missing genes 2 > 1" in mapped.quality.reasons


def test_render_locus_map_json_is_a_separate_document() -> None:
    """Given a detailed report, When rendered, Then it is gspec.locus_map/1
    with the file + map payload and no compact cluster/2 fields."""
    alignment = build_locus_alignment(LOCUS, [_row(":500")])
    mapped = build_locus_map(_call(), alignment, FEATURES, (_call(),))
    report = ClusterReport(
        file="sample.fna", best=None, loci=(_call(),), locus_maps={"OL1": mapped}
    )
    rendered = render_locus_map_json(
        [report], ClusterParams(db="vpa_ok"), now=datetime(2026, 11, 2, tzinfo=UTC)
    )
    assert '"schema": "gspec.locus_map/1"' in rendered
    assert '"file": "sample.fna"' in rendered
    assert '"locus": "OL1"' in rendered
    assert '"quality"' in rendered


def test_extra_panel_genes_require_family_novelty_and_query_overlap() -> None:
    """Given another locus gene whose family is absent from the best locus
    and whose span overlaps the best alignment, When mapped, Then it is
    reported as an extra panel gene."""
    other_locus = LocusFeatures(
        id="OL2",
        label="O2",
        type="O2",
        genes=(
            GeneFeature(gene_id="b1", start=1, end=100, strand="+"),
            GeneFeature(gene_id="b2", start=201, end=300, strand="+"),
        ),
    )
    features = FeaturesDocument(loci=(LOCUS, other_locus))
    families = GeneFamilyDocument(
        families=(
            GeneFamily(family_id="GF1", genes=("g1", "g2", "g3", "b1")),
            GeneFamily(family_id="GF2", genes=("b2",)),
        ),
        loci=(
            LocusFamilySet(locus="OL1", families=("GF1",)),
            LocusFamilySet(locus="OL2", families=("GF1", "GF2")),
        ),
        containment=(),
    )
    other_call = LocusCall(
        locus="OL2",
        label="O2",
        type="O2",
        coverage_pct=99.0,
        identity_pct=99.0,
        rank=2,
        genes=(
            GeneCall(
                gene_id="b1",
                start=1,
                end=100,
                strand="+",
                coverage_pct=100.0,
                identity_pct=100.0,
                verdict="present",
            ),
            GeneCall(
                gene_id="b2",
                start=201,
                end=300,
                strand="+",
                coverage_pct=100.0,
                identity_pct=100.0,
                verdict="present",
            ),
        ),
        missing=(),
    )
    alignment = build_locus_alignment(LOCUS, [_row(":500")])
    mapped = build_locus_map(
        _call(),
        alignment,
        features,
        (_call(), other_call),
        families,
        {"OL2": {"b2": (QuerySpan(contig="ctg1", qstart=100, qend=200),)}},
    )
    assert [
        (gene.gene_id, gene.family, gene.source_locus) for gene in mapped.extra_panel_genes
    ] == [("b2", "GF2", "OL2")]
    assert mapped.quality.extra_panel_genes == 1


def test_screen_cluster_file_detail_attaches_locus_map(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the engine's PAF rows, When detail is requested, Then the report
    carries the rich map while the default path stays map-free."""
    database = Database(
        name="test",
        path=tmp_path,
        sequences_path=tmp_path / "sequences",
        kind="cluster",
    )
    row = PafRecord(
        qname="ctg1",
        qlen=600,
        qstart=0,
        qend=500,
        strand="+",
        tname="OL1",
        tlen=500,
        tstart=0,
        tend=500,
        nmatch=500,
        alen=500,
        mapq=60,
        cs=":500",
    )

    def _fake_rows(*args: object, **kwargs: object) -> Iterator[PafRecord]:
        return iter([row])

    monkeypatch.setattr("gspec.cluster.iter_minimap2", _fake_rows)
    params = ClusterParams(db="test")
    detailed = screen_cluster_file(tmp_path / "sample.fna", database, FEATURES, params, detail=True)
    assert "OL1" in detailed.locus_maps
    assert detailed.locus_maps["OL1"].quality.status == "pass"
    assert detailed.locus_maps["OL1"].insertions == ()
    plain = screen_cluster_file(tmp_path / "sample.fna", database, FEATURES, params)
    assert plain.locus_maps == {}
