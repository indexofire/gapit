"""Unit tests for the abricate-byte-compatible TSV/CSV formatter (SPEC.md §5)."""

from gapit.blast import BlastRow
from gapit.formats.tsv import HEADER, format_tsv
from gapit.hits import process_rows
from gapit.report import Report

TET_ID = "tinyamr~~~tetA~~~NC_000913.3:100-900~~~TETRACYCLINE"


def row(
    qseqid: str = "contig1",
    qstart: int = 1,
    qend: int = 79,
    qlen: int = 79,
    sseqid: str = TET_ID,
    sstart: int = 1,
    send: int = 79,
    slen: int = 79,
    sstrand: str = "plus",
    evalue: float = 1e-40,
    length: int = 79,
    pident: float = 100.0,
    gaps: int = 0,
    gapopen: int = 0,
    stitle: str = TET_ID + " tetracycline efflux pump TetA",
) -> BlastRow:
    """A full-length 100% tetA hit with per-test overrides."""
    return BlastRow(
        qseqid=qseqid,
        qstart=qstart,
        qend=qend,
        qlen=qlen,
        sseqid=sseqid,
        sstart=sstart,
        send=send,
        slen=slen,
        sstrand=sstrand,
        evalue=evalue,
        length=length,
        pident=pident,
        gaps=gaps,
        gapopen=gapopen,
        stitle=stitle,
    )


def make_report(rows: list[BlastRow], file: str = "in.fa") -> Report:
    """Process rows like screen_file does (mincov=0 keeps everything)."""
    hits = process_rows(rows, mincov=0.0, default_db="tinyamr")
    return Report(file=file, hits=tuple(sorted(hits, key=lambda hit: (hit.sequence, hit.start))))


FULL_LINE = (
    "in.fa\tcontig1\t1\t79\t+\ttetA\t1-79/79\t===============\t0/0\t"
    "100.00\t100.00\ttinyamr\tNC_000913.3:100-900\ttetracycline efflux pump TetA\tTETRACYCLINE"
)


def test_header_present_by_default() -> None:
    """Given one hit, When formatted, Then the output is HEADER line + one row."""
    output = format_tsv([make_report([row()])], csv=False, noheader=False, nopath=False)
    assert output == HEADER + "\n" + FULL_LINE + "\n"


def test_header_absent_with_noheader() -> None:
    """Given noheader, When formatted, Then only hit lines, no '#FILE' line."""
    output = format_tsv([make_report([row()])], csv=False, noheader=True, nopath=False)
    assert output == FULL_LINE + "\n"
    assert "#FILE" not in output


def test_empty_reports_noheader_is_empty_string() -> None:
    """Given no reports and noheader, When formatted, Then the output is ''."""
    assert format_tsv([], csv=False, noheader=True, nopath=False) == ""


def test_empty_hits_still_print_header() -> None:
    """Given reports with zero hits, When formatted, Then the header alone."""
    assert format_tsv([make_report([])], csv=False, noheader=False, nopath=False) == HEADER + "\n"


def test_csv_separator_everywhere() -> None:
    """Given csv=True, When formatted, Then header and rows use ',' — and no
    field ever contains a ',' (product cleanup strips them), so join is safe."""
    output = format_tsv([make_report([row()])], csv=True, noheader=False, nopath=False)
    assert output == HEADER.replace("\t", ",") + "\n" + FULL_LINE.replace("\t", ",") + "\n"


def test_nopath_basenames_the_file_column() -> None:
    """Given nopath and a deep file path, When formatted, Then FILE is basename."""
    report = make_report([row()], file="/deep/dir/in.fa")
    output = format_tsv([report], csv=False, noheader=True, nopath=True)
    assert output.startswith("in.fa\tcontig1\t")


def test_coverage_and_gaps_string_shapes() -> None:
    """Given s_start/s_end/s_len 1-94/94 and gap_openings/gaps 1/3, When
    formatted, Then COVERAGE is '1-94/94' and GAPS is '1/3'."""
    hits = process_rows(
        [row(slen=94, length=97, gaps=3, gapopen=1, send=94, qend=97, qlen=97)],
        mincov=0.0,
        default_db="tinyamr",
    )
    report = Report(file="x", hits=tuple(hits))
    output = format_tsv([report], csv=False, noheader=True, nopath=True)
    assert "\t1-94/94\t" in output
    assert "\t1/3\t" in output


def test_coverage_79_996_renders_as_80_00() -> None:
    """Given coverage_pct=79.996 (SPEC §4: displays 80.00 even though it is
    filtered at mincov=80; here kept at mincov=0), When formatted, Then
    %COVERAGE shows '80.00'."""
    hits = process_rows(
        [row(length=20000, gaps=1, slen=25000, qend=20000, qlen=20000)],
        mincov=0.0,
        default_db="tinyamr",
    )
    report = Report(file="x", hits=tuple(hits))
    output = format_tsv([report], csv=False, noheader=True, nopath=True)
    assert "\t80.00\t100.00\t" in output


def test_identity_rounds_to_two_decimals() -> None:
    """Given pident 96.907, When formatted, Then %IDENTITY is '96.91'."""
    hits = process_rows([row(pident=96.907)], mincov=0.0, default_db="tinyamr")
    report = Report(file="x", hits=tuple(hits))
    output = format_tsv([report], csv=False, noheader=True, nopath=True)
    assert "\t100.00\t96.91\t" in output


def test_multi_report_single_header_and_order_preserved() -> None:
    """Given two reports, When formatted, Then one header and report order
    preserved (format does not re-sort; screen_file sorted within reports)."""
    first = make_report([row()], file="b.fa")
    second = make_report([row(qseqid="contig0")], file="a.fa")
    output = format_tsv([first, second], csv=False, noheader=False, nopath=False)
    lines = output.splitlines()
    assert len(lines) == 3
    assert lines[0] == HEADER
    assert lines[1].startswith("b.fa\tcontig1")
    assert lines[2].startswith("a.fa\tcontig0")
    assert output.endswith("\n") and not output.endswith("\n\n")
