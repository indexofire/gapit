"""Cluster-build boundary: locus/gene fields must not break table output."""

from pathlib import Path

import pytest

from gapit.clusterbuild import build_cluster_database
from gapit.errors import DatabaseError
from gapit.gbfeatures import GeneFeature, LocusFeatures, LocusSequence


def _locus(gene_id: str) -> LocusSequence:
    return LocusSequence(
        features=LocusFeatures(
            id="L1",
            label="L1",
            type="L1",
            genes=(GeneFeature(gene_id=gene_id, start=1, end=4, strand="+"),),
        ),
        sequence="ACGT",
    )


def test_cluster_build_rejects_control_char_in_gene_id(tmp_path: Path) -> None:
    """Given a locus whose gene id contains a tab, When built, Then
    BUILD_INVALID fires before any artifact is written."""
    db_dir = tmp_path / "cluster"
    db_dir.mkdir()

    with pytest.raises(DatabaseError) as excinfo:
        build_cluster_database(
            db_dir,
            name="cluster",
            loci=[_locus("g\t1")],
            typing=None,
            source_urls=(),
            fetched_at="2026-10-22T00:00:00Z",
        )

    assert excinfo.value.code == "BUILD_INVALID"
    assert excinfo.value.context["locus"] == "L1"
    assert excinfo.value.context["field"] == "gene_id"
    assert not (db_dir / "sequences").exists()
    assert not (db_dir / "features.json").exists()
