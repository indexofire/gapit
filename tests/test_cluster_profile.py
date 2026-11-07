"""Per-database cluster profile tests: opt-in coverage, thresholds, fallback
evidence, and locus-status semantics."""

from __future__ import annotations

from pathlib import Path

import pytest

from gspec.bundled import find_bundled
from gspec.cluster import runs_by_locus
from gspec.cluster_profile import (
    CLUSTER_PROFILE_FILENAME,
    load_cluster_profile,
)
from gspec.db import Database
from gspec.errors import DatabaseError
from gspec.paf import PafRecord
from gspec.screening_cluster import resolve_cluster_params


def _database(tmp_path: Path) -> Database:
    return Database(
        name="cluster-test",
        path=tmp_path,
        sequences_path=tmp_path / "sequences",
        kind="cluster",
    )


def _row(*, primary: bool) -> PafRecord:
    return PafRecord(
        qname="contig",
        qlen=100,
        qstart=0,
        qend=10,
        strand="+",
        tname="OL1",
        tlen=10,
        tstart=0,
        tend=10,
        nmatch=10,
        alen=10,
        mapq=60,
        is_primary=primary,
        cs=":10",
    )


def test_missing_profile_keeps_historical_defaults(tmp_path: Path) -> None:
    """Given no cluster.json, When the profile is loaded, Then primary-only
    and the historical thresholds are in effect."""
    profile = load_cluster_profile(tmp_path)
    assert profile.coverage_mode == "primary"
    assert profile.min_cluster_cov is None
    params = resolve_cluster_params(_database(tmp_path), 90, 90, 96, 1)
    assert params.coverage_mode == "primary"
    assert params.min_cluster_cov == 96


def test_profile_opts_into_secondary_and_thresholds(tmp_path: Path) -> None:
    """Given a cluster profile, When params are resolved, Then the database
    opts into secondary coverage and profile defaults, while a non-default
    CLI threshold still wins."""
    (tmp_path / CLUSTER_PROFILE_FILENAME).write_text(
        '{"schema":"gspec.cluster_profile/1","coverage_mode":"primary_secondary",'
        '"min_gene_cov":50.0,"min_cluster_cov":90.0,"fallback_evidence":true,'
        '"locus_typeable_coverage":50.0,"locus_typeable_identity":90.0,'
        '"locus_quality_min_coverage":80.0,"locus_quality_min_identity":90.0,'
        '"locus_quality_max_missing_genes":3,"locus_quality_max_extra_panel_genes":2}',
        encoding="utf-8",
    )
    db = _database(tmp_path)
    params = resolve_cluster_params(db, 90, 90, 96, 1)
    assert params.coverage_mode == "primary_secondary"
    assert params.fallback_evidence is True
    assert params.min_gene_cov == 50.0
    assert params.min_cluster_cov == 90.0
    assert params.locus_typeable_coverage == 50.0
    assert params.locus_typeable_identity == 90.0
    assert params.locus_quality_min_coverage == 80.0
    assert params.locus_quality_min_identity == 90.0
    assert params.locus_quality_max_missing_genes == 3
    assert params.locus_quality_max_extra_panel_genes == 2
    explicit = resolve_cluster_params(db, 60, 90, 95, 1)
    assert explicit.min_gene_cov == 60.0
    assert explicit.min_cluster_cov == 95.0


def test_malformed_profile_is_typed(tmp_path: Path) -> None:
    """Given malformed cluster.json, When loaded, Then fail with a typed
    database error rather than silently changing coverage semantics."""
    (tmp_path / CLUSTER_PROFILE_FILENAME).write_text("{not json", encoding="utf-8")
    with pytest.raises(DatabaseError) as raised:
        load_cluster_profile(tmp_path)
    assert raised.value.code == "CLUSTER_PROFILE_MALFORMED"


def test_unknown_profile_key_is_typed(tmp_path: Path) -> None:
    """Given a misspelled profile key, When loaded, Then fail instead of
    silently ignoring a policy switch."""
    (tmp_path / CLUSTER_PROFILE_FILENAME).write_text(
        '{"schema":"gspec.cluster_profile/1","fallback_evidance":true}',
        encoding="utf-8",
    )
    with pytest.raises(DatabaseError) as raised:
        load_cluster_profile(tmp_path)
    assert raised.value.code == "CLUSTER_PROFILE_MALFORMED"


def test_secondary_mode_includes_secondary_rows(tmp_path: Path) -> None:
    """Given a secondary-only alignment, When cluster math is primary-only,
    Then it is ignored; the profile's primary_secondary mode includes it."""
    rows = [_row(primary=False)]
    assert runs_by_locus(rows) == {}
    grouped = runs_by_locus(rows, coverage_mode="primary_secondary")
    assert "OL1" in grouped
    assert grouped["OL1"][1] == 10


def test_bundled_vpa_ok_opts_into_secondary_coverage() -> None:
    """Given the bundled vpa_ok snapshot, When its profile is loaded, Then
    it asks for Kaptive-style primary+secondary coverage."""
    bundled = find_bundled("vpa_ok")
    assert bundled is not None and bundled.profile_path is not None
    profile = load_cluster_profile(bundled.path)
    assert profile.coverage_mode == "primary_secondary"
    assert profile.fallback_evidence is True
    assert profile.min_gene_cov == 50.0
    assert profile.locus_typeable_coverage == 50.0
    assert profile.locus_typeable_identity == 90.0
    assert profile.locus_quality_min_coverage == 80.0
    assert profile.locus_quality_min_identity == 90.0
    assert profile.locus_quality_max_missing_genes == 3
    assert profile.locus_quality_max_extra_panel_genes == 2


def test_fallback_evidence_is_scheme_family_scoped() -> None:
    """Given a low-confidence O/K fallback, When evidence is enabled, Then
    O reports the best OL locus and K reports the best KL locus — never a
    cross-family locus."""
    from gspec.cluster import ClusterReport, LocusCall
    from gspec.typing_cluster2 import evaluate_cluster2
    from gspec.typing_models import TypingDocument

    ol = LocusCall(
        locus="OL3",
        label="OL3",
        type="OL3",
        coverage_pct=99.0,
        identity_pct=99.0,
        rank=1,
        genes=(),
        missing=(),
    )
    kl = LocusCall(
        locus="KL29",
        label="KL29",
        type="KL29",
        coverage_pct=30.0,
        identity_pct=98.0,
        rank=2,
        genes=(),
        missing=(),
    )
    report = ClusterReport(file="sample.fna", best=None, loci=(ol, kl))
    rule: dict[str, object] = {
        "model": "cluster_match",
        "phenotype": "O1",
        "coverage": {"weight": 0.55, "floor": 90.0},
        "identity": {"weight": 0.25, "floor": 95.0},
        "key_genes": {"weight": 0.2, "genes": []},
    }
    document = TypingDocument.model_validate(
        {
            "schemes": [
                {
                    "name": "o_group",
                    "rules": [rule],
                    "cutoff": 0.9,
                    "ambiguity_margin": 0.05,
                    "fallback": "OUT",
                },
                {
                    "name": "k_group",
                    "rules": [{**rule, "phenotype": "K1"}],
                    "cutoff": 0.9,
                    "ambiguity_margin": 0.05,
                    "fallback": "KUT",
                },
            ]
        }
    )
    calls = evaluate_cluster2(
        report,
        document,
        fallback_evidence=True,
        typeable_coverage=50.0,
        typeable_identity=90.0,
    )
    assert calls["o_group"].winner_locus == "OL3"
    assert calls["k_group"].winner_locus == "KL29"
    assert calls["o_group"].locus_status == "unknown_locus"
    assert calls["k_group"].locus_status == "unknown_locus"
    assert calls["o_group"].locus_confidence == "typeable"
    assert calls["k_group"].locus_confidence == "untypeable"
    plain = evaluate_cluster2(report, document)
    assert plain["o_group"].winner_locus is None
    assert plain["k_group"].winner_locus is None
    assert plain["o_group"].locus_status is None
    assert plain["o_group"].locus_confidence is None
