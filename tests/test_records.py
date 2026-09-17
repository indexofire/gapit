"""Tests for the truth-source Record JSONL store and the gapit.manifest/1 sidecar."""

from pathlib import Path

import pytest
from pydantic import ValidationError

from gapit.errors import InputError
from gapit.formats.json import ToolDocument
from gapit.records import (
    Manifest,
    Record,
    count_records,
    read_manifest,
    read_records,
    write_manifest,
    write_records,
)

DIGEST = "a3f5d0" + "0" * 58
RECORDS: tuple[Record, ...] = (
    Record(
        db="tinyamr",
        gene="tetA",
        sequence="ATGCGTTTAC",
        accession="A00001",
        function=("tetracycline", "doxycycline"),
        product="β-lactamase-like efflux pump",
        source_id="NUH0001",
    ),
    Record(db="tinyamr", gene="blaSHV", sequence="ATGGCA"),
)


def make_manifest(*, sha256: str = DIGEST, n_records: int = 2) -> Manifest:
    return Manifest(
        name="tinyamr",
        source_urls=("https://example.org/tinyamr.fasta", "https://mirror.example/tinyamr"),
        fetched_at="2026-09-17T10:30:00Z",
        sha256=sha256,
        n_records=n_records,
        dbtype="nucl",
    )


def test_jsonl_round_trip_preserves_values_including_defaults(tmp_path: Path) -> None:
    """Given one fully populated record and one all-defaults record, When
    written then read back, Then every field value survives exactly."""
    path = tmp_path / "records.jsonl"
    write_records(RECORDS, path)
    assert list(read_records(path)) == list(RECORDS)
    minimal = RECORDS[1]
    assert minimal.accession == ""
    assert minimal.function == ()
    assert minimal.product == "n/a"
    assert minimal.source_id == ""


def test_write_is_deterministic(tmp_path: Path) -> None:
    """Given the same records in the same order, When written to two files,
    Then the bytes are identical: one JSON object per line, LF endings only."""
    first = tmp_path / "a.jsonl"
    second = tmp_path / "b.jsonl"
    write_records(RECORDS, first)
    write_records(RECORDS, second)
    payload = first.read_bytes()
    assert payload == second.read_bytes()
    assert b"\r" not in payload
    assert payload.endswith(b"\n")
    assert len(payload.splitlines()) == len(RECORDS)


def test_write_empty_iterable_makes_empty_file(tmp_path: Path) -> None:
    """Given no records, When written, Then the file is empty and reads back
    as zero records."""
    path = tmp_path / "records.jsonl"
    write_records((), path)
    assert path.read_bytes() == b""
    assert list(read_records(path)) == []


def test_streaming_read_yields_early_records_before_error(tmp_path: Path) -> None:
    """Given a file whose first record is valid and a later line is garbage,
    When streamed, Then the valid record is yielded first and the failure
    carries the physical line number in its context."""
    path = tmp_path / "records.jsonl"
    path.write_text(RECORDS[0].model_dump_json() + "\n\n{not json\n", encoding="utf-8")
    stream = read_records(path)
    assert next(stream) == RECORDS[0]
    with pytest.raises(InputError) as excinfo:
        next(stream)
    assert excinfo.value.code == "RECORDS_MALFORMED"
    assert excinfo.value.context == {"file": str(path), "line": "3"}


def test_malformed_json_line_raises_typed_error_with_line_number(tmp_path: Path) -> None:
    """Given a line that is not valid JSON, When the file is read, Then an
    InputError RECORDS_MALFORMED is raised with file/line context and the
    pydantic reason chained and embedded in the message."""
    path = tmp_path / "records.jsonl"
    path.write_text(RECORDS[0].model_dump_json() + "\n{not json\n", encoding="utf-8")
    with pytest.raises(InputError) as excinfo:
        list(read_records(path))
    error = excinfo.value
    assert error.code == "RECORDS_MALFORMED"
    assert error.context == {"file": str(path), "line": "2"}
    assert str(error).startswith(f"malformed record at {path}:2:")
    assert isinstance(error.__cause__, ValidationError)
    assert " ".join(str(error.__cause__).split()) in str(error)


def test_schema_violating_line_raises_malformed(tmp_path: Path) -> None:
    """Given a line of valid JSON missing a required field, When read, Then
    the same RECORDS_MALFORMED error fires with that line number."""
    path = tmp_path / "records.jsonl"
    path.write_text('{"db": "tinyamr"}\n', encoding="utf-8")
    with pytest.raises(InputError) as excinfo:
        list(read_records(path))
    assert excinfo.value.code == "RECORDS_MALFORMED"
    assert excinfo.value.context == {"file": str(path), "line": "1"}


def test_missing_records_file_raises_input_not_found(tmp_path: Path) -> None:
    """Given a nonexistent records file, When read, Then InputError
    INPUT_NOT_FOUND mirrors the cmd_db wording style with the path."""
    missing = tmp_path / "nope.jsonl"
    with pytest.raises(InputError) as excinfo:
        list(read_records(missing))
    error = excinfo.value
    assert error.code == "INPUT_NOT_FOUND"
    assert error.context == {"file": str(missing)}
    assert "not found or unreadable" in str(error)


def test_blank_lines_are_skipped(tmp_path: Path) -> None:
    """Given blank lines around and between records, When read, Then exactly
    the records come back in order."""
    path = tmp_path / "records.jsonl"
    path.write_text(
        "\n" + RECORDS[0].model_dump_json() + "\n\n" + RECORDS[1].model_dump_json() + "\n\n",
        encoding="utf-8",
    )
    assert list(read_records(path)) == list(RECORDS)


def test_count_records_matches_with_blank_lines(tmp_path: Path) -> None:
    """Given a records file with interleaved blank lines, When counted, Then
    only record lines are counted."""
    path = tmp_path / "records.jsonl"
    path.write_text(
        "\n".join(["", RECORDS[0].model_dump_json(), "", RECORDS[1].model_dump_json(), ""]),
        encoding="utf-8",
    )
    assert count_records(path) == len(RECORDS)


def test_manifest_round_trip_with_schema_alias_first(tmp_path: Path) -> None:
    """Given a manifest, When written then read, Then the bytes are indented
    JSON with the `schema` alias as first key and a trailing newline, and the
    model round-trips equal."""
    path = tmp_path / "gapit-manifest.json"
    manifest = make_manifest()
    write_manifest(manifest, path)
    text = path.read_text(encoding="utf-8")
    assert text.startswith('{\n  "schema": "gapit.manifest/1",')
    assert text.endswith("}\n")
    assert read_manifest(path) == manifest


def test_manifest_rejects_short_sha256() -> None:
    """Given a sha256 that is not 64 hex chars, When constructing, Then
    pydantic validation fails."""
    with pytest.raises(ValidationError):
        make_manifest(sha256="deadbeef")


def test_manifest_normalizes_uppercase_sha256_to_lowercase() -> None:
    """Given an uppercase 64-hex digest, When constructed, Then it is stored
    lowercase."""
    assert make_manifest(sha256=DIGEST.upper()).sha256 == DIGEST


def test_manifest_rejects_negative_n_records() -> None:
    """Given n_records < 0, When constructing, Then validation fails."""
    with pytest.raises(ValidationError):
        make_manifest(n_records=-1)


def test_manifest_defaults_applied() -> None:
    """Given a manifest constructed without the optional fields, Then the
    frozen-contract defaults hold."""
    manifest = make_manifest()
    assert manifest.schema_name == "gapit.manifest/1"
    assert manifest.header_format == "gapit/v1"
    assert manifest.upstream_version == ""
    assert manifest.tool == ToolDocument()
    assert manifest.makeblastdb_version == ""
    assert manifest.minimap2_version == ""


def test_read_manifest_missing_file_raises_input_not_found(tmp_path: Path) -> None:
    missing = tmp_path / "nope.json"
    with pytest.raises(InputError) as excinfo:
        read_manifest(missing)
    assert excinfo.value.code == "INPUT_NOT_FOUND"
    assert excinfo.value.context == {"file": str(missing)}


def test_read_manifest_malformed_raises_typed_error(tmp_path: Path) -> None:
    """Given a manifest with a wrong-typed field, When read, Then InputError
    MANIFEST_MALFORMED carries the file context."""
    path = tmp_path / "gapit-manifest.json"
    path.write_text('{"schema": "gapit.manifest/1", "name": 42}', encoding="utf-8")
    with pytest.raises(InputError) as excinfo:
        read_manifest(path)
    error = excinfo.value
    assert error.code == "MANIFEST_MALFORMED"
    assert error.context == {"file": str(path)}
