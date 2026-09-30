"""Tests for the GFF3 cluster parser (gapit.gffparse).

Embedded ``##FASTA`` sections and sidecar ``<stem>.fa`` files must both work;
each FASTA sequence becomes one locus and CDS rows become its genes. GFF
coordinates are 1-based inclusive and pass through unchanged.
"""

from pathlib import Path

import pytest

from gapit.errors import InputError
from gapit.gffparse import parse_gff3_features

DATA = Path(__file__).parent / "data" / "cluster"

SCL1 = "CATTGCAAGGTTACGATCGG" * 9  # 180 bp
SCL2 = "GGTACCTTAGGCATCGAATT" * 6  # 120 bp


def test_embedded_fasta_yields_loci_and_genes() -> None:
    """Given a GFF3 with an embedded ##FASTA section, When parsed, Then each
    FASTA sequence is one locus and its CDS rows are the genes — coordinates
    kept 1-based inclusive, strands carried through, non-CDS rows ignored."""
    loci = parse_gff3_features(DATA / "clusters.gff3")

    assert [locus.features.id for locus in loci] == ["cl1", "cl2"]
    cl1, cl2 = (locus.features for locus in loci)
    assert (cl1.label, cl1.type) == ("cl1", "cl1")
    wzx, wzy, man_c = cl1.genes
    assert (wzx.gene_id, wzx.start, wzx.end, wzx.strand) == ("wzxK", 1, 60, "+")
    assert wzx.product == "export protein"
    assert (wzy.gene_id, wzy.start, wzy.end, wzy.strand) == ("cl1_002", 61, 120, "-")
    assert (man_c.gene_id, man_c.start, man_c.end, man_c.strand) == ("manC", 121, 180, "+")

    (big_a,) = cl2.genes
    assert (big_a.gene_id, big_a.start, big_a.end, big_a.strand) == ("bigA", 10, 90, "+")
    assert big_a.product == "big protein"


def test_embedded_fasta_sequences_extracted() -> None:
    """Given the embedded-FASTA fixture, When parsed, Then each locus carries
    the full sequence of its FASTA record."""
    loci = parse_gff3_features(DATA / "clusters.gff3")

    assert loci[0].sequence == SCL1
    assert loci[1].sequence == SCL2


def test_sidecar_fasta_resolves_and_parses() -> None:
    """Given a GFF3 without ##FASTA but with a <stem>.fa sidecar, When
    parsed, Then the sidecar supplies the locus sequences."""
    loci = parse_gff3_features(DATA / "clusters_sidecar.gff3")

    assert [locus.features.id for locus in loci] == ["cl1", "cl2"]
    assert [len(locus.features.genes) for locus in loci] == [3, 1]
    assert loci[0].sequence == SCL1


def test_gzipped_gff_with_embedded_fasta_parses() -> None:
    """Given a .gz GFF3, When parsed, Then decompression is transparent."""
    loci = parse_gff3_features(DATA / "clusters_gz.gff3.gz")

    assert [locus.features.id for locus in loci] == ["cl1", "cl2"]


def test_missing_fasta_raises_listing_searched_paths() -> None:
    """Given a GFF3 with neither ##FASTA nor any sidecar, When parsed, Then a
    typed InputError lists every sidecar path that was searched."""
    with pytest.raises(InputError) as excinfo:
        parse_gff3_features(DATA / "nofasta.gff3")

    assert excinfo.value.code == "GFF_WITHOUT_FASTA"
    searched = excinfo.value.context["searched"].split(";")
    assert any(path.endswith("nofasta.fa") for path in searched)
    assert any(path.endswith("nofasta.fasta") for path in searched)


def test_wrong_gff_version_is_rejected() -> None:
    """Given a file whose first directive is not ##gff-version 3, When
    parsed, Then a typed InputError names the offending line."""
    bad = DATA / "clusters.gff3"
    text = bad.read_text(encoding="utf-8").replace("##gff-version 3", "##gff-version 2", 1)
    stray = bad.with_name("stray.gff3")
    stray.write_text(text, encoding="utf-8")

    try:
        with pytest.raises(InputError) as excinfo:
            parse_gff3_features(stray)
        assert excinfo.value.code == "INVALID_GFF"
        assert "##gff-version 3" in str(excinfo.value)
    finally:
        stray.unlink()


def test_cds_on_unknown_sequence_is_rejected() -> None:
    """Given a GFF3 whose CDS row names a sequence the FASTA does not carry,
    When parsed, Then a typed InputError names the orphan sequence."""
    text = (DATA / "clusters.gff3").read_text(encoding="utf-8")
    head, _, fasta = text.partition("##FASTA\n")
    stray = DATA / "orphan.gff3"
    stray.write_text(
        f"{head}ghost\tbakta\tCDS\t1\t50\t.\t+\t0\tID=g1;gene=boo\n##FASTA\n{fasta}",
        encoding="utf-8",
    )

    try:
        with pytest.raises(InputError) as excinfo:
            parse_gff3_features(stray)
        assert excinfo.value.code == "INVALID_GFF"
        assert excinfo.value.context["seqid"] == "ghost"
    finally:
        stray.unlink()
