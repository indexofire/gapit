"""CSV rendering must be RFC 4180-compatible, not comma-joined TSV.

These tests are adversarial on purpose: gene/accession/function/file/sample
values are user- or database-controlled, and product cleanup only removes
commas from the product field. A CSV parser must see the same number of
fields as the header, with the original values intact.
"""

import csv
import io
from pathlib import Path

from gapit.formats.cluster import format_cluster_tsv
from gapit.formats.reads_tsv import format_reads_tsv
from gapit.formats.summary import format_summary_tsv
from gapit.formats.tsv import format_tsv
from gapit.hits import Hit
from gapit.reads import GeneCoverage, ReadsReport
from gapit.report import Report
from gapit.summary import SummaryMatrix, SummaryParams, SummaryRow


def _hit(**overrides: object) -> Hit:
    values: dict[str, object] = {
        "sequence": "ctg",
        "start": 1,
        "end": 9,
        "strand": "+",
        "gene": "a,b",
        "database": "db",
        "accession": 'A"1',
        "product": "p,q",
        "function": "class,one",
        "s_start": 1,
        "s_end": 9,
        "s_len": 9,
        "coverage_map": "===============",
        "gap_openings": 0,
        "gaps": 0,
        "identity_pct": 99.0,
        "coverage_pct": 100.0,
        "aligned_len": 9,
    }
    values.update(overrides)
    return Hit.model_validate(values)


def _rows(text: str) -> list[list[str]]:
    return list(csv.reader(io.StringIO(text)))


def test_contig_csv_quotes_commas_and_double_quotes() -> None:
    """Given fields containing commas/quotes, When rendered as CSV, Then every
    row has 15 parsed fields and the values round-trip exactly."""
    report = Report(file="/tmp/in,put.fa", hits=(_hit(),))
    rows = _rows(format_tsv([report], csv=True, noheader=False, nopath=False))

    assert len(rows[0]) == 15
    assert len(rows[1]) == 15
    assert rows[1][0] == "/tmp/in,put.fa"
    assert rows[1][5] == "a,b"
    assert rows[1][12] == 'A"1'
    assert rows[1][13] == "p,q"
    assert rows[1][14] == "class,one"


def test_paired_reads_csv_keeps_the_sample_field_whole() -> None:
    """Given a paired sample key ``s1, s2``, When rendered as CSV, Then the
    header and row both have nine fields and the sample stays one cell."""
    gene = GeneCoverage(
        database="db",
        gene="g,1",
        accession="a",
        function="f",
        product="p",
        tlen=100,
        breadth_pct=99.0,
        mean_depth=3.0,
        reads_mapped=5,
        present=True,
        mean_identity_pct=99.0,
    )
    report = ReadsReport(reads=("s1", "s2"), genes=(gene,))
    rows = _rows(format_reads_tsv([report], csv=True))

    assert len(rows[0]) == 9
    assert len(rows[1]) == 9
    assert rows[1][0] == "s1, s2"
    assert rows[1][1] == "g,1"


def test_summary_csv_quotes_a_comma_gene_in_the_header() -> None:
    """Given a gene name containing a comma, When rendered as a summary CSV,
    Then the header and row both have three fields and the gene is intact."""
    matrix = SummaryMatrix(
        params=SummaryParams(),
        genes=("a,b",),
        rows=(SummaryRow(file="in,fa", num_found=1, cells={"a,b": ("100.00",)}),),
    )
    rows = _rows(format_summary_tsv(matrix, csv=True))

    assert len(rows[0]) == 3
    assert len(rows[1]) == 3
    assert rows[0][2] == "a,b"
    assert rows[1][0] == "in,fa"


def test_cluster_csv_keeps_file_and_row_widths_in_sync() -> None:
    """Given a file path containing a comma, When rendered as cluster CSV,
    Then the header and row parse to the same field count."""
    from gapit.cluster import ClusterReport

    report = ClusterReport(file="/tmp/in,put.fa", best=None, loci=())
    rows = _rows(format_cluster_tsv([report], csv=True, noheader=False, nopath=False, typed=False))

    assert len(rows[0]) == len(rows[1])
    assert rows[1][0] == "/tmp/in,put.fa"


def test_cli_golden_csv_file_still_uses_lf(tmp_path: Path) -> None:
    """Given a simple CSV render, When written, Then the line terminator is
    LF (not the csv module's default CRLF) so golden files stay stable."""
    report = Report(
        file="in.fa", hits=(_hit(gene="g1", accession="A1", product="p", function="f"),)
    )
    text = format_tsv([report], csv=True, noheader=False, nopath=False)
    assert "\r" not in text
    assert text.endswith("\n")


def test_summary_reads_its_own_quoted_csv() -> None:
    """Given a CSV screen table with a comma-bearing file/gene, When summary
    parses stdin, Then the quoted values survive as single cells."""
    from gapit.summary import build_summary

    text = (
        "#FILE,SEQUENCE,START,END,STRAND,GENE,COVERAGE,COVERAGE_MAP,GAPS,"
        "%COVERAGE,%IDENTITY,DATABASE,ACCESSION,PRODUCT,RESISTANCE\n"
        '"in,fa",ctg,1,9,+,"a,b",1-9/9,===============,0/0,100.00,99.00,db,A,p,f\n'
    )
    matrix = build_summary([], SummaryParams(), warn=lambda message: None, stdin_text=text)

    assert matrix.genes == ("a,b",)
    assert matrix.rows[0].file == "in,fa"


def test_typing_input_reads_its_own_quoted_csv() -> None:
    """Given a quoted CSV screen table, When typing parses stdin, Then the
    quoted FILE/GENE cells are not split."""
    from gapit.typing_input import read_typing_input

    text = (
        "#FILE,SEQUENCE,START,END,STRAND,GENE,COVERAGE,COVERAGE_MAP,GAPS,"
        "%COVERAGE,%IDENTITY,DATABASE,ACCESSION,PRODUCT,RESISTANCE\n"
        '"in,fa",ctg,1,9,+,"a,b",1-9/9,===============,0/0,100.00,99.00,db,A,p,f\n'
    )
    parsed = read_typing_input([], stdin_text=text)

    assert parsed.files[0].file == "in,fa"
    assert "a,b" in parsed.files[0].calls
