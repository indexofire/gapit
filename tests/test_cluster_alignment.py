"""Query/target cs walk + LocusAlignment evidence tests."""

from __future__ import annotations

from typing import Literal

from gspec.cluster_alignment import build_locus_alignment, gene_query_spans
from gspec.cluster_math import parse_cs, parse_cs_steps
from gspec.gbfeatures import GeneFeature, LocusFeatures
from gspec.paf import PafRecord

LOCUS = LocusFeatures(
    id="OL1",
    label="O1",
    type="O1",
    genes=(
        GeneFeature(gene_id="g1", start=1, end=100, strand="+"),
        GeneFeature(gene_id="g2", start=201, end=300, strand="+"),
        GeneFeature(gene_id="g3", start=401, end=500, strand="+"),
    ),
)


def _row(
    *,
    strand: Literal["+", "-"] = "+",
    cs: str,
    nmatch: int = 750,
    alen: int = 750,
    primary: bool = True,
) -> PafRecord:
    return PafRecord(
        qname="ctg1",
        qlen=1000,
        qstart=100,
        qend=850,
        strand=strand,
        tname="OL1",
        tlen=500,
        tstart=0,
        tend=500,
        nmatch=nmatch,
        alen=alen,
        mapq=60,
        is_primary=primary,
        cs=cs,
    )


def test_parse_cs_steps_keeps_query_and_target_walks() -> None:
    """Given insertions/deletions, When parsed, Then both coordinate spaces
    advance correctly while the legacy target runs stay unchanged."""
    cs = f":100+{'A' * 150}:100-{'C' * 50}:250"
    steps = parse_cs_steps(cs)
    assert [(s.op, s.qstart, s.qend, s.tstart, s.tend) for s in steps] == [
        ("match", 0, 100, 0, 100),
        ("insertion", 100, 250, 100, 100),
        ("match", 250, 350, 100, 200),
        ("deletion", 350, 350, 200, 250),
        ("match", 350, 600, 250, 500),
    ]
    assert parse_cs(cs, 0) == ((0, 100, 100), (100, 200, 100), (200, 250, 0), (250, 500, 250))


def test_forward_insertions_report_contig_interval_and_flanks() -> None:
    """Given two large insertions, When mapped, Then each carries its
    forward query interval, target boundary, and flanking locus genes."""
    alignment = build_locus_alignment(
        LOCUS,
        [_row(cs=f":100+{'A' * 150}:200+{'C' * 100}:200")],
    )
    assert alignment.length == 500
    assert len(alignment.pieces) == 1
    assert [(i.qstart, i.qend, i.length, i.target_position) for i in alignment.insertions] == [
        (200, 350, 150, 100),
        (550, 650, 100, 300),
    ]
    assert [(i.left_gene, i.right_gene, i.flank) for i in alignment.insertions] == [
        ("g1", "g2", "between"),
        ("g2", "g3", "between"),
    ]


def test_reverse_insertions_use_forward_contig_coordinates() -> None:
    """Given a reverse PAF row, When an insertion is mapped, Then its query
    interval is transformed back into forward contig coordinates."""
    alignment = build_locus_alignment(
        LOCUS,
        [_row(strand="-", cs=f":100+{'A' * 150}:200+{'C' * 100}:200")],
    )
    assert [(i.qstart, i.qend) for i in alignment.insertions] == [(600, 750), (300, 400)]


def test_insertion_inside_gene_is_marked() -> None:
    """Given an insertion inside a reference gene, When flanked, Then the
    within-gene signal is explicit instead of forcing a between-genes call."""
    alignment = build_locus_alignment(LOCUS, [_row(cs=f":50+{'A' * 120}:450")])
    (insertion,) = alignment.insertions
    assert insertion.within_gene == "g1"
    assert insertion.flank == "within_gene"


def test_deletion_reports_spanned_genes() -> None:
    """Given a deletion covering g2, When mapped, Then the deletion names
    the spanned gene and its flanking genes."""
    alignment = build_locus_alignment(LOCUS, [_row(cs=f":100-{'A' * 150}:250")])
    (deletion,) = alignment.deletions
    assert (deletion.tstart, deletion.tend, deletion.length) == (100, 250, 150)
    assert deletion.spanned_genes == ("g2",)
    assert (deletion.left_gene, deletion.right_gene) == ("g1", "g3")


def test_small_indels_are_omitted_but_pieces_remain() -> None:
    """Given a below-threshold indel, When mapped, Then the piece stays and
    the indel is omitted."""
    alignment = build_locus_alignment(LOCUS, [_row(cs=f":100+{'A' * 10}:400")], min_indel_bp=100)
    assert len(alignment.pieces) == 1
    assert alignment.insertions == ()


def test_primary_mode_excludes_secondary_rows() -> None:
    """Given a secondary-only row, When primary mode is requested, Then no
    evidence is emitted; primary_secondary includes it."""
    row = _row(cs=":500", primary=False)
    assert build_locus_alignment(LOCUS, [row]).pieces == ()
    assert len(build_locus_alignment(LOCUS, [row], coverage_mode="primary_secondary").pieces) == 1


def test_gene_query_spans_project_only_aligned_reference_bases() -> None:
    """Given insertions between genes, When gene spans are projected, Then
    each span tracks the reference gene's query interval, not the insertion."""
    spans = gene_query_spans(LOCUS, [_row(cs=f":100+{'A' * 150}:200+{'C' * 100}:200")])
    assert [(s.contig, s.qstart, s.qend) for s in spans["g1"]] == [("ctg1", 100, 200)]
    assert [(s.qstart, s.qend) for s in spans["g2"]] == [(450, 550)]
    assert [(s.qstart, s.qend) for s in spans["g3"]] == [(750, 850)]
