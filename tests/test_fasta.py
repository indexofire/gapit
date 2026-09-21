"""Unit tests for the streaming FASTA reader (gapit.fasta)."""

import bz2
import gzip
from pathlib import Path

import pytest
from pydantic import ValidationError

from gapit.errors import InputError
from gapit.fasta import FastaRecord, iter_fasta, iter_fasta_headers

FIXTURE = Path(__file__).parent / "data" / "db" / "tinyamr" / "sequences"


def write_fasta(path: Path, content: str) -> Path:
    """Write plain FASTA content to path."""
    path.write_text(content, encoding="utf-8")
    return path


def test_iter_fasta_joins_multiline_records(tmp_path: Path) -> None:
    """Given a two-record multi-line FASTA, When parsed, Then sequences are the
    concatenated lines and id/description split on the first whitespace."""
    path = write_fasta(tmp_path / "two.fa", ">a desc one\nACGT\nTGCA\n>b\nGGCC\n")
    records = list(iter_fasta(path))
    assert records == [
        FastaRecord(id="a", description="desc one", sequence="ACGTTGCA"),
        FastaRecord(id="b", description="", sequence="GGCC"),
    ]


def test_iter_fasta_preserves_sequence_case(tmp_path: Path) -> None:
    """Given lowercase sequence letters, When parsed, Then case is preserved as-is."""
    records = list(iter_fasta(write_fasta(tmp_path / "lc.fa", ">x\nacgTn\n")))
    assert records[0].sequence == "acgTn"


def test_iter_fasta_reads_the_committed_fixture() -> None:
    """Given the tinyamr fixture, When parsed, Then 3 records with canonical headers
    and ACGT-only sequences come out."""
    records = list(iter_fasta(FIXTURE))
    assert [r.id for r in records] == [
        "tinyamr~~~tetA~~~NC_000913.3:100-900~~~TETRACYCLINE",
        "tinyamr~~~blaTEM-1~~~J01749.1:1-861~~~BETA-LACTAM",
        "tinyamr~~~sul1~~~U12338.4:1-940~~~SULFONAMIDE",
    ]
    assert records[0].description == "tetracycline efflux pump TetA"
    assert all(set(r.sequence) <= set("ACGT") for r in records)


def test_iter_fasta_reads_gz(tmp_path: Path) -> None:
    """Given the fixture gzip-compressed, When parsed, Then records equal the plain ones."""
    content = FIXTURE.read_text(encoding="utf-8")
    gz_path = tmp_path / "db.fa.gz"
    gz_path.write_bytes(gzip.compress(content.encode("utf-8")))
    assert list(iter_fasta(gz_path)) == list(iter_fasta(FIXTURE))


def test_iter_fasta_reads_bz2(tmp_path: Path) -> None:
    """Given the fixture bzip2-compressed, When parsed, Then records equal the plain ones."""
    content = FIXTURE.read_text(encoding="utf-8")
    bz2_path = tmp_path / "db.fa.bz2"
    bz2_path.write_bytes(bz2.compress(content.encode("utf-8")))
    assert list(iter_fasta(bz2_path)) == list(iter_fasta(FIXTURE))


def test_iter_fasta_rejects_content_before_first_header(tmp_path: Path) -> None:
    """Given junk before the first '>', When parsed, Then InputError (exit 5)."""
    path = write_fasta(tmp_path / "junk.fa", "ACGT\n>a\nACGT\n")
    with pytest.raises(InputError) as excinfo:
        list(iter_fasta(path))
    assert excinfo.value.exit_code == 5


def test_iter_fasta_rejects_record_with_empty_sequence(tmp_path: Path) -> None:
    """Given a header immediately followed by another header (empty sequence),
    When parsed, Then InputError."""
    path = write_fasta(tmp_path / "empty.fa", ">a\nACGT\n>b\n")
    with pytest.raises(InputError):
        list(iter_fasta(path))


def test_iter_fasta_empty_file_yields_no_records(tmp_path: Path) -> None:
    """Given an empty file, When parsed, Then zero records and no error."""
    assert list(iter_fasta(write_fasta(tmp_path / "empty.fa", ""))) == []


def test_fasta_record_is_frozen(tmp_path: Path) -> None:
    """Given a parsed record, When a field is assigned, Then pydantic rejects it
    (setattr because the typed model already marks fields read-only)."""
    record = next(iter(iter_fasta(write_fasta(tmp_path / "one.fa", ">a\nACGT\n"))))
    attribute = "id"  # variable keeps ruff B010 quiet; the model marks fields read-only
    with pytest.raises(ValidationError):
        setattr(record, attribute, "mutated")


def test_iter_fasta_headers_extracts_id_and_description(tmp_path: Path) -> None:
    """Given a multi-record FASTA, When headers are iterated, Then (id,
    description) pairs come out with the description split on the first
    whitespace — same split as iter_fasta."""
    path = write_fasta(tmp_path / "two.fa", ">a desc one\nACGT\nTGCA\n>b\nGGCC\n")
    assert list(iter_fasta_headers(path)) == [("a", "desc one"), ("b", "")]


def test_iter_fasta_headers_matches_iter_fasta_on_fixture() -> None:
    """Given the committed tinyamr fixture, When both iterators run, Then the
    header-only pairs equal iter_fasta's (id, description) pairs exactly."""
    assert list(iter_fasta_headers(FIXTURE)) == [(r.id, r.description) for r in iter_fasta(FIXTURE)]


def test_iter_fasta_headers_skips_sequence_lines(tmp_path: Path) -> None:
    """Given records with long multi-line sequences, When headers are
    iterated, Then sequence lines are skipped and never surface."""
    path = write_fasta(tmp_path / "seqs.fa", ">a d\n" + "ACGT\n" * 1000 + ">b e\nGGCC\n")
    assert list(iter_fasta_headers(path)) == [("a", "d"), ("b", "e")]


def test_iter_fasta_headers_allows_empty_sequence(tmp_path: Path) -> None:
    """Given a header immediately followed by another header, When headers
    are iterated, Then NO InputError — the empty-sequence check does not
    apply because sequences are never read."""
    path = write_fasta(tmp_path / "empty.fa", ">a\nACGT\n>b\n")
    assert list(iter_fasta_headers(path)) == [("a", ""), ("b", "")]


def test_iter_fasta_headers_rejects_content_before_first_header(tmp_path: Path) -> None:
    """Given junk before the first '>', When headers are iterated, Then the
    same InputError as iter_fasta (malformed FASTA is still rejected)."""
    path = write_fasta(tmp_path / "junk.fa", "ACGT\n>a\nACGT\n")
    with pytest.raises(InputError) as excinfo:
        list(iter_fasta_headers(path))
    assert excinfo.value.code == "INVALID_FASTA"


def test_iter_fasta_headers_empty_file_yields_nothing(tmp_path: Path) -> None:
    """Given an empty file, When headers are iterated, Then zero pairs."""
    assert list(iter_fasta_headers(write_fasta(tmp_path / "empty.fa", ""))) == []


def test_iter_fasta_headers_reads_gz(tmp_path: Path) -> None:
    """Given a gzip-wrapped FASTA, When headers are iterated, Then pairs
    equal the plain file's (shared transparent decompression)."""
    plain = write_fasta(tmp_path / "plain.fa", ">a desc one\nACGT\n>b\nGGCC\n")
    gz_path = tmp_path / "headers.fa.gz"
    gz_path.write_bytes(gzip.compress(plain.read_bytes()))
    assert list(iter_fasta_headers(gz_path)) == [("a", "desc one"), ("b", "")]


def test_iter_fasta_truncated_gz_raises_invalid_fasta(tmp_path: Path) -> None:
    """Given a .gz FASTA cut to ~60% of its bytes, When iterated, Then the
    mid-iteration EOFError becomes InputError INVALID_FASTA with the file in
    context (typed, not an uncaught EOFError escaping the generator)."""
    blob = gzip.compress(b">a desc\n" + b"ACGT" * 100 + b"\n>b\nGGCC\n")
    gz_path = tmp_path / "truncated.fa.gz"
    gz_path.write_bytes(blob[: len(blob) * 3 // 5])
    with pytest.raises(InputError) as excinfo:
        list(iter_fasta(gz_path))
    assert excinfo.value.code == "INVALID_FASTA"
    assert excinfo.value.context == {"file": str(gz_path)}


def test_iter_fasta_headers_truncated_gz_raises_invalid_fasta(tmp_path: Path) -> None:
    """Given the same truncated .gz, When headers are iterated, Then the same
    typed InputError INVALID_FASTA fires (both readers share the wrap)."""
    blob = gzip.compress(b">a desc\n" + b"ACGT" * 100 + b"\n>b\nGGCC\n")
    gz_path = tmp_path / "truncated.fa.gz"
    gz_path.write_bytes(blob[: len(blob) * 3 // 5])
    with pytest.raises(InputError) as excinfo:
        list(iter_fasta_headers(gz_path))
    assert excinfo.value.code == "INVALID_FASTA"
    assert excinfo.value.context == {"file": str(gz_path)}
