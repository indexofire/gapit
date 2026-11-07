"""Gene-family sidecar and containment tests (pure; no BLAST)."""

from __future__ import annotations

from pathlib import Path

import pytest

from gspec.errors import DatabaseError
from gspec.gbfeatures import FeaturesDocument, GeneFeature, LocusFeatures
from gspec.gene_families import (
    GeneFamilyDocument,
    HomologyHit,
    build_family_document,
    cluster_families,
    filter_homology_hits,
    read_gene_families,
    write_gene_families,
)

FEATURES = FeaturesDocument(
    loci=(
        LocusFeatures(
            id="A",
            label="A",
            type="A",
            genes=(
                GeneFeature(gene_id="A1", start=1, end=100, strand="+"),
                GeneFeature(gene_id="A2", start=201, end=300, strand="+"),
            ),
        ),
        LocusFeatures(
            id="B",
            label="B",
            type="B",
            genes=(
                GeneFeature(gene_id="B1", start=1, end=100, strand="+"),
                GeneFeature(gene_id="B2", start=201, end=300, strand="+"),
                GeneFeature(gene_id="B3", start=401, end=500, strand="+"),
            ),
        ),
        LocusFeatures(
            id="C",
            label="C",
            type="C",
            genes=(GeneFeature(gene_id="C1", start=1, end=100, strand="+"),),
        ),
    )
)


def _hit(query: str, subject: str, *, identity: float = 95.0, alen: int = 100) -> HomologyHit:
    return HomologyHit(
        query=query,
        subject=subject,
        identity_pct=identity,
        alignment_len=alen,
        query_len=100,
        subject_len=100,
        evalue=1e-20,
    )


def test_filter_homology_hits_uses_identity_and_shorter_coverage() -> None:
    """Given a low-identity and a short-coverage hit, When filtered, Then
    only the family-level hit survives."""
    hits = [
        _hit("A1", "B1"),
        _hit("A2", "B2", identity=70.0),
        _hit("A1", "C1", alen=50),
    ]
    kept = filter_homology_hits(hits)
    assert [(h.query, h.subject) for h in kept] == [("A1", "B1")]


def test_cluster_families_is_deterministic_and_merges_components() -> None:
    """Given A1-B1 and B1-B2 edges, When clustered, Then all three genes
    share one family and ids follow the lexicographic first member."""
    assignments = cluster_families(
        ["B2", "A1", "B1", "C1"],
        [_hit("A1", "B1"), _hit("B1", "B2")],
    )
    assert assignments["A1"] == assignments["B1"] == assignments["B2"]
    assert assignments["C1"] != assignments["A1"]
    assert assignments["A1"] == "GF000001"


def test_build_document_records_strict_containment_and_diagnostics() -> None:
    """Given A's families are a strict subset of B's, When assembled, Then
    the containment edge carries B's diagnostic family."""
    assignments = {
        "A1": "GF000001",
        "A2": "GF000002",
        "B1": "GF000001",
        "B2": "GF000002",
        "B3": "GF000003",
        "C1": "GF000004",
    }
    document = build_family_document(FEATURES, assignments)
    edges = {(edge.subset, edge.superset): edge.diagnostic for edge in document.containment}
    assert edges == {("A", "B"): ("GF000003",)}
    assert document.families[0].genes == ("A1", "B1")


def test_build_document_rejects_missing_assignment() -> None:
    """Given a partially assigned panel, When assembled, Then fail instead
    of silently treating the missing gene as unique."""
    with pytest.raises(DatabaseError) as raised:
        build_family_document(FEATURES, {"A1": "GF000001"})
    assert raised.value.code == "GENE_FAMILIES_INCOMPLETE"


def test_round_trip_is_byte_stable(tmp_path: Path) -> None:
    """Given a family document, When written and read, Then the model
    round-trips and the serialized file is stable."""
    assignments = {
        "A1": "GF000001",
        "A2": "GF000002",
        "B1": "GF000001",
        "B2": "GF000002",
        "B3": "GF000003",
        "C1": "GF000004",
    }
    document = build_family_document(FEATURES, assignments)
    path = tmp_path / "gene_families.json"
    write_gene_families(document, path)
    first = path.read_bytes()
    write_gene_families(read_gene_families(path), path)
    assert path.read_bytes() == first
    assert isinstance(read_gene_families(path), GeneFamilyDocument)
