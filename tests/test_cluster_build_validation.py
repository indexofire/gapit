"""Cluster-build boundary: locus/gene fields must not break table output."""

from pathlib import Path

import pytest

from gapit.clusterbuild import build_cluster_database
from gapit.errors import DatabaseError
from gapit.gbfeatures import GeneFeature, LocusFeatures, LocusSequence


def _locus(gene_id: str, *, locus_id: str = "L1") -> LocusSequence:
    return LocusSequence(
        features=LocusFeatures(
            id=locus_id,
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


def test_cluster_build_rejects_duplicate_locus_ids(tmp_path: Path) -> None:
    """Given two loci with the same id, When built, Then BUILD_INVALID fires
    before any artifact is written (the engine keys on the locus id)."""
    db_dir = tmp_path / "cluster"
    db_dir.mkdir()

    with pytest.raises(DatabaseError) as excinfo:
        build_cluster_database(
            db_dir,
            name="cluster",
            loci=[_locus("g1"), _locus("g2")],
            typing=None,
            source_urls=(),
            fetched_at="2026-10-22T00:00:00Z",
        )

    assert excinfo.value.code == "BUILD_INVALID"
    assert excinfo.value.context["reason"] == "duplicate"
    assert excinfo.value.context["locus"] == "L1"
    assert not (db_dir / "sequences").exists()


def test_cluster_build_rejects_table_delimiter_in_name(tmp_path: Path) -> None:
    """Given a cluster db name with a tab, When built, Then BUILD_INVALID
    fires before any artifact is written."""
    db_dir = tmp_path / "cluster"
    db_dir.mkdir()

    with pytest.raises(DatabaseError) as excinfo:
        build_cluster_database(
            db_dir,
            name="bad\tname",
            loci=[_locus("g1")],
            typing=None,
            source_urls=(),
            fetched_at="2026-10-22T00:00:00Z",
        )

    assert excinfo.value.code == "BUILD_INVALID"
    assert excinfo.value.context["db"] == "bad\tname"
    assert not (db_dir / "sequences").exists()


def test_cluster_build_rejects_duplicate_gene_ids_within_a_locus(tmp_path: Path) -> None:
    """Given two CDS sharing one gene id inside a locus, When built, Then
    BUILD_INVALID fires before any artifact is written."""
    locus = LocusSequence(
        features=LocusFeatures(
            id="L1",
            label="L1",
            type="L1",
            genes=(
                GeneFeature(gene_id="dup", start=1, end=4, strand="+"),
                GeneFeature(gene_id="dup", start=5, end=8, strand="+"),
            ),
        ),
        sequence="ACGTACGT",
    )
    db_dir = tmp_path / "cluster"
    db_dir.mkdir()

    with pytest.raises(DatabaseError) as excinfo:
        build_cluster_database(
            db_dir,
            name="cluster",
            loci=[locus],
            typing=None,
            source_urls=(),
            fetched_at="2026-10-22T00:00:00Z",
        )

    assert excinfo.value.code == "BUILD_INVALID"
    assert excinfo.value.context["reason"] == "duplicate"
    assert not (db_dir / "sequences").exists()
