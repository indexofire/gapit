"""Unit tests for the streaming FASTA reader (gaita.fasta)."""

import bz2
import gzip
from pathlib import Path

import pytest
from pydantic import ValidationError

from gaita.errors import InputError
from gaita.fasta import FastaRecord, iter_fasta

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
