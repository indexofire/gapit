"""Tests for the GenBank FEATURES/ORIGIN parser (gapit.gbfeatures).

Fixtures under tests/data/cluster cover both real-world styles the cluster
databases ship in: Bakta/modern (LOCUS id + CDS qualifiers) and kaptive-style
(source ``/note`` label pairs). Malformed inputs must raise typed InputErrors.
"""

from pathlib import Path

import pytest

from gapit.errors import InputError
from gapit.gbfeatures import parse_genbank_features

DATA = Path(__file__).parent / "data" / "cluster"
REAL_CPS = Path("/tmp/opencode/vpautils/src/vpautils/database/serotype/CPSgc.gbk")

S9001 = ("ACGTTGCAAGGCTTACGATC" * 9).upper()  # 180 bp
S9002 = ("TTGACGATCGGTAACCAGTT" * 10).upper()  # 200 bp


def test_bakta_style_loci_and_genes() -> None:
    """Given a Bakta-style two-locus GBK, When parsed, Then two loci come
    back with the LOCUS-name ids, the label from the source ``K locus`` note,
    the type falling back to the locus id (no type note), and per-locus gene
    lists of three genes each."""
    loci = parse_genbank_features(DATA / "bakta_style.gbk")

    assert [locus.features.id for locus in loci] == ["KL9001", "KL9002"]
    first, second = (locus.features for locus in loci)
    assert first.label == "KL9001"
    assert first.type == "KL9001"
    assert second.label == "KL9002"
    assert second.type == "KL9002"
    assert len(first.genes) == 3
    assert len(second.genes) == 3


def test_bakta_style_gene_coordinates_strand_and_names() -> None:
    """Given the Bakta-style fixture, When parsed, Then CDS coordinates are
    1-based inclusive with +/- strands, the gene name prefers /gene over
    /locus_tag, and products carry through."""
    first = parse_genbank_features(DATA / "bakta_style.gbk")[0].features
    wzx, wzy, man_c = first.genes

    assert (wzx.gene_id, wzx.start, wzx.end, wzx.strand) == ("wzx", 1, 60, "-")
    assert wzx.product == "capsule polysaccharide export protein"
    assert (wzy.gene_id, wzy.start, wzy.end, wzy.strand) == ("KL9001_002", 61, 120, "+")
    assert wzy.product == "capsule polysaccharide polymerase"
    assert (man_c.gene_id, man_c.start, man_c.end, man_c.strand) == ("manC", 121, 180, "-")

    second = parse_genbank_features(DATA / "bakta_style.gbk")[1].features
    orf_a, orf_b, orf_c = second.genes
    assert (orf_a.gene_id, orf_a.start, orf_a.end, orf_a.strand) == ("orfA", 5, 70, "+")
    assert (orf_b.gene_id, orf_b.start, orf_b.end, orf_b.strand) == ("orfB", 80, 140, "-")
    assert (orf_c.gene_id, orf_c.start, orf_c.end, orf_c.strand) == ("KL9002_009", 150, 200, "+")


def test_bakta_style_extracts_uppercased_origin_sequences() -> None:
    """Given lowercase ORIGIN data, When parsed, Then each locus carries its
    whole uppercased sequence (the seqconvert normalization contract)."""
    loci = parse_genbank_features(DATA / "bakta_style.gbk")

    assert loci[0].sequence == S9001
    assert loci[1].sequence == S9002


def test_kaptive_style_label_and_type_from_source_notes() -> None:
    """Given a kaptive-style GBK, When parsed, Then label and type come from
    the ``K locus``/``K type`` source notes — including the literal
    ``unknown`` type."""
    loci = parse_genbank_features(DATA / "kaptive_style.gbk")

    assert [(locus.features.label, locus.features.type) for locus in loci] == [
        ("KL106", "K6"),
        ("KL107", "unknown"),
    ]


def test_unnamed_cds_fall_back_to_note_then_positional_id() -> None:
    """Given a kaptive-style record whose CDS carry /gene, /note-only orf
    names, or no name at all (the real O-locus database shape), When parsed,
    Then gene ids resolve /gene -> /note -> the positional <locus>_NN id —
    no gene is dropped for being unnamed."""
    (locus,) = parse_genbank_features(DATA / "unnamed_cds.gbk")

    assert [gene.gene_id for gene in locus.features.genes] == [
        "wzm",
        "orf11",
        "AB819963_03",
    ]


def test_gzipped_input_parses() -> None:
    """Given a .gz GBK, When parsed, Then the same structured records come
    back through the transparent decompressor, with ``<``/````>``
    partial-location prefixes stripped."""
    loci = parse_genbank_features(DATA / "gz_input.gbk.gz")

    assert len(loci) == 1
    genes = loci[0].features.genes
    assert [(g.gene_id, g.start, g.end, g.strand) for g in genes] == [
        ("gzA", 1, 60, "+"),
        ("gzB", 61, 120, "-"),
    ]


def test_compound_join_cds_is_rejected() -> None:
    """Given a CDS with a compound join() location, When parsed, Then a typed
    InputError names the locus and the unsupported location (unsupported
    compound CDS — documented limitation)."""
    with pytest.raises(InputError) as excinfo:
        parse_genbank_features(DATA / "compound.gbk")

    assert excinfo.value.code == "UNSUPPORTED_LOCATION"
    assert excinfo.value.context["locus"] == "KLJOIN"
    assert "join(" in excinfo.value.context["location"]


@pytest.mark.skipif(not REAL_CPS.is_file(), reason="real CPSgc.gbk not present (offline CI)")
def test_real_bakta_database_parses_163_loci() -> None:
    """Given the rightsholder's real 16 MB Bakta-style CPSgc.gbk, When
    parsed, Then exactly 163 loci come back, every locus carries genes, and
    every gene sits inside its locus bounds."""
    loci = parse_genbank_features(REAL_CPS)

    assert len(loci) == 163
    assert all(len(locus.features.genes) > 0 for locus in loci)
    for locus in loci:
        length = len(locus.sequence)
        assert length > 0
        for gene in locus.features.genes:
            assert 1 <= gene.start <= gene.end <= length, (locus.features.id, gene.gene_id)
