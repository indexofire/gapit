"""Unit tests for gapit.seqconvert (native any2fasta -q -u replacement).

Semantics locked here mirror the perl (see the seqconvert.py docstring for
line references); the differential suite (test_seqconvert_differential.py,
parity env on PATH) proves parsed-record equality against the real binary.
"""

import bz2
import gzip
from pathlib import Path

import pytest

from gapit.errors import InputError
from gapit.seqconvert import SeqFormat, detect_format, to_fasta_lines

CONVERT = Path(__file__).parent / "data" / "convert"
TETASEQ = "ATGGTCAATTCCGCTGGCGCTGGTTATCGGTGGCACTGGCTGCGTTTTGGCGATGGTTTTCGGCAACCGTGCGCTGGCA"


def convert(path: Path) -> str:
    return "".join(to_fasta_lines(path))


def parsed(text: str) -> list[tuple[str, str, str]]:
    """(id, description, sequence) triples from FASTA text."""
    records: list[tuple[str, str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append(_split(header, "".join(chunks)))
            header, chunks = line[1:], []
        else:
            chunks.append(line)
    if header is not None:
        records.append(_split(header, "".join(chunks)))
    return records


def _split(header: str, sequence: str) -> tuple[str, str, str]:
    parts = header.split(None, 1)
    return (parts[0] if parts else "", parts[1] if len(parts) > 1 else "", sequence)


def test_fasta_passthrough_uppercases_and_keeps_headers_verbatim() -> None:
    """Given the committed mixed-case FASTA fixture (with a blank line), When
    converted, Then headers are verbatim (id + description) and the sequence
    is uppercased with blank lines skipped."""
    records = parsed(convert(CONVERT / "sample.fa"))
    assert records == [
        ("FAKESEQ01", "synthetic fixture sequence one", TETASEQ + "GTTAGGCCAAT"),
        ("FAKESEQ02", "", "TTTAAACCC"),
    ]


def test_fasta_wraps_sequence_at_60_columns(tmp_path: Path) -> None:
    """Given a 130-nt lowercase FASTA record, When converted, Then sequence
    lines are 60/60/10 uppercase columns."""
    query = tmp_path / "long.fa"
    query.write_text(">long_one\n" + "ta" * 65 + "\n", encoding="utf-8")
    lines = convert(query).splitlines()
    assert lines[0] == ">long_one"
    assert [len(line) for line in lines[1:]] == [60, 60, 10]
    assert all(line == "TA" * (len(line) // 2) for line in lines[1:])


def test_fastq_uses_full_header_line_and_uppercases() -> None:
    """Given the committed FASTQ fixture, When converted, Then each 4-line
    record becomes a FASTA record: '@'-stripped header kept verbatim (id +
    description), sequence uppercased, '+' and quality dropped."""
    records = parsed(convert(CONVERT / "sample.fq"))
    assert records == [
        ("FAKEREAD01", "synthetic read number one", "ACGTNACGTNACGTNACGTN"),
        ("FAKEREAD02", "synthetic read number two", "GGGGCCCCAAAATTTT"),
    ]


def test_fastq_trailing_lone_header_emits_nothing(tmp_path: Path) -> None:
    """Given a FASTQ whose last record is an unterminated '@' header, When
    converted, Then no record is emitted for it (perl loop guard i < $#lines)."""
    query = tmp_path / "trailing.fq"
    query.write_text("@r1 desc\nacgt\n+\nIIII\n@r9 alone\n", encoding="utf-8")
    records = parsed(convert(query))
    assert records == [("r1", "desc", "ACGT")]


def test_genbank_multi_record_locus_id_and_uppercase() -> None:
    """Given the committed multi-record GBK fixture, When converted, Then ids
    are the LOCUS names, DEFINITION contributes no description (perl never
    reads it), VERSION is ignored (no -g), and ORIGIN coordinates are
    stripped while the sequence is uppercased."""
    records = parsed(convert(CONVERT / "sample.gbk"))
    assert records == [
        ("FAKECTG01", "", TETASEQ),
        ("FAKECTG02", "", "C" * 60),
    ]


def test_genbank_keeps_digits_past_column_10(tmp_path: Path) -> None:
    """Given an ORIGIN data line whose sequence field contains digits, When
    converted, Then the digits survive (perl strips only whitespace after the
    10-column coordinate prefix, unlike the EMBL parser)."""
    query = tmp_path / "digits.gbk"
    query.write_text(
        "LOCUS       DGTS999 6 bp    DNA     linear   UNK 01-JAN-2026\n"
        "ORIGIN\n"
        "        1 ab12cd\n"
        "//\n",
        encoding="utf-8",
    )
    assert parsed(convert(query)) == [("DGTS999", "", "AB12CD")]


def test_genbank_unterminated_record_is_dropped(tmp_path: Path) -> None:
    """Given a GBK record never closed by //, When converted, Then it is not
    emitted (perl flushes only at //)."""
    query = tmp_path / "unterminated.gbk"
    query.write_text(
        "LOCUS       LOST999 4 bp    DNA     linear   UNK 01-JAN-2026\nORIGIN\n        1 acgt\n",
        encoding="utf-8",
    )
    assert convert(query) == ""


def test_genbank_empty_record_emits_header_only(tmp_path: Path) -> None:
    """Given ORIGIN immediately followed by //, When converted, Then the
    record is a bare header line (perl prints an empty sequence)."""
    query = tmp_path / "empty.gbk"
    query.write_text(
        "LOCUS       VOID999 0 bp    DNA     linear   UNK 01-JAN-2026\nORIGIN\n//\n",
        encoding="utf-8",
    )
    assert convert(query) == ">VOID999\n"


def test_embl_multi_record_id_up_to_semicolon() -> None:
    """Given the committed multi-record EMBL fixture, When converted, Then ids
    are the ID-line text up to the first ';', DE contributes no description,
    and SQ digits/whitespace are stripped with the sequence uppercased."""
    records = parsed(convert(CONVERT / "sample.embl"))
    assert records == [
        ("FAKECTG01", "", TETASEQ),
        ("FAKECTG02", "", "C" * 60),
    ]


def test_gzipped_input_matches_plain(tmp_path: Path) -> None:
    """Given the committed .gz fixtures, When converted, Then output equals
    the plain-file conversion (transparent decompression)."""
    assert convert(CONVERT / "sample.fa.gz") == convert(CONVERT / "sample.fa")
    assert convert(CONVERT / "sample.gbk.gz") == convert(CONVERT / "sample.gbk")


def test_bzip2_input_matches_plain(tmp_path: Path) -> None:
    """Given the committed .bz2 fixture, When converted, Then output equals
    the plain-file conversion."""
    assert convert(CONVERT / "sample.embl.bz2") == convert(CONVERT / "sample.embl")


def test_gz_written_in_test_also_works(tmp_path: Path) -> None:
    """Given a freshly gzipped FASTA, When converted, Then it parses (the
    opener keys off the .gz suffix)."""
    query = tmp_path / "made.fa.gz"
    query.write_bytes(gzip.compress(b">mk1\nacgt\n"))
    assert parsed(convert(query)) == [("mk1", "", "ACGT")]


def test_detect_format_matrix(tmp_path: Path) -> None:
    """Given first lines for each format, When detected, Then the format
    matches the perl's regex table (LOCUS/ID need a space; >/@ need a
    non-whitespace character after the marker)."""
    cases = [
        (">seq1 some description\nACGT\n", SeqFormat.fasta),
        ("@read1 desc\nACGT\n+\nII\n", SeqFormat.fastq),
        ("LOCUS       abc 79 bp\nORIGIN\n//\n", SeqFormat.genbank),
        ("ID   abc; SV 1;\nSQ   Sequence 4 BP;\n//\n", SeqFormat.embl),
    ]
    for text, expected in cases:
        probe = tmp_path / f"probe{expected.value}"
        probe.write_text(text, encoding="utf-8")
        assert detect_format(probe) is expected


def test_detect_format_rejects_marker_then_space(tmp_path: Path) -> None:
    """>Given '> ' and '@ ' first lines (marker + whitespace), When detected,
    Then InputError fires: perl requires ^>\\S, a header must start with a
    non-whitespace character."""
    for marker in (">", "@"):
        probe = tmp_path / f"space{marker.strip()}"
        probe.write_text(f"{marker} leading space\n", encoding="utf-8")
        with pytest.raises(InputError):
            detect_format(probe)


def test_empty_file_raises_invalid_input(tmp_path: Path) -> None:
    """Given an empty file, When converted, Then InputError INVALID_INPUT
    with the perl's empty-input message."""
    query = tmp_path / "empty.fa"
    query.write_text("", encoding="utf-8")
    with pytest.raises(InputError) as excinfo:
        convert(query)
    assert excinfo.value.code == "INVALID_INPUT"
    assert "appears to be empty" in str(excinfo.value)


def test_junk_first_line_raises_invalid_input(tmp_path: Path) -> None:
    """Given a non-sequence first line, When converted, Then InputError
    INVALID_INPUT carrying the perl's (typo'd) unfamiliar-format message."""
    query = tmp_path / "junk.txt"
    query.write_text("this is not sequence data at all\n", encoding="utf-8")
    with pytest.raises(InputError) as excinfo:
        convert(query)
    assert excinfo.value.code == "INVALID_INPUT"
    assert str(excinfo.value).startswith("Unfamilar format with first line:")


def test_corrupt_gzip_raises_invalid_input(tmp_path: Path) -> None:
    """Given a .gz file that is not gzip at all, When converted, Then the
    decompression failure becomes InputError INVALID_INPUT (typed, not an
    uncaught OSError)."""
    query = tmp_path / "junk.fa.gz"
    query.write_bytes(b"this is not gzip data")
    with pytest.raises(InputError) as excinfo:
        convert(query)
    assert excinfo.value.code == "INVALID_INPUT"


def test_truncated_gz_raises_invalid_input(tmp_path: Path) -> None:
    """Given a valid .gz cut to ~60% of its bytes (partial download), When
    converted, Then the mid-read EOFError — not an OSError subclass — becomes
    InputError INVALID_INPUT (detect_format passes on the first readable line;
    the failure surfaces mid-iteration)."""
    blob = gzip.compress(b">mk1\n" + b"acgt" * 100 + b"\n>mk2\n" + b"tttt" * 100 + b"\n")
    query = tmp_path / "truncated.fa.gz"
    query.write_bytes(blob[: len(blob) * 3 // 5])
    with pytest.raises(InputError) as excinfo:
        convert(query)
    assert excinfo.value.code == "INVALID_INPUT"


def test_truncated_bz2_raises_invalid_input(tmp_path: Path) -> None:
    """Given a valid .bz2 cut to ~60% of its bytes, When converted, Then the
    same InputError INVALID_INPUT — bz2 also raises EOFError on a truncated
    stream, not an OSError."""
    blob = bz2.compress(b">mk1\n" + b"acgt" * 100 + b"\n>mk2\n" + b"tttt" * 100 + b"\n")
    query = tmp_path / "truncated.fa.bz2"
    query.write_bytes(blob[: len(blob) * 3 // 5])
    with pytest.raises(InputError) as excinfo:
        convert(query)
    assert excinfo.value.code == "INVALID_INPUT"
