"""Tests for the Wave A3 build pipeline: records.jsonl -> built database.

Binary-requiring tests run the pixi env's real makeblastdb/blastn/minimap2
offline (same pattern as test_screen_integration.py)."""

import hashlib
from pathlib import Path

import pytest

from gapit import __version__, dbbuild
from gapit.dbbuild import build_database, generate_sequences, verify_sequences
from gapit.errors import DatabaseError
from gapit.records import Record, read_manifest, write_records

DB = "tinyamr"
SEQ_WRAPPED = "A" * 70  # exercises the 60-col wrap: 60 + 10
SEQ_SHORT = "ACGTACGTACGTACGTACGT"
SEQ_TAIL = "AC" * 30 + "G"  # 61 chars: 60 + 1

RECORDS: tuple[Record, ...] = (
    Record(
        db=DB,
        gene="ant(3'')-Ia|v2",
        sequence=SEQ_WRAPPED,
        accession="A00001",
        function=("streptomycin", "spectinomycin"),
        product="aminoglycoside adenylyltransferase",
        source_id="NUH0001",
    ),
    Record(
        db=DB,
        gene="mcr-1",
        sequence=SEQ_SHORT,
        accession="MCR00002",
        function=("colistin", "polymyxin"),
        product="phosphoethanolamine transferase",
    ),
    Record(db=DB, gene="tetA", sequence=SEQ_TAIL),
)

# Hand-computed gapit/v1 headers: values percent-escaped by the dbcodec rules
# (| -> %7C, ; -> %3B), empty fields render as empty segments, defaults apply.
HEADER_WRAPPED = (
    ">gapit|db=tinyamr|gene=ant(3'')-Ia%7Cv2|acc=A00001"
    "|func=streptomycin%3Bspectinomycin aminoglycoside adenylyltransferase"
)
HEADER_SHORT = (
    ">gapit|db=tinyamr|gene=mcr-1|acc=MCR00002"
    "|func=colistin%3Bpolymyxin phosphoethanolamine transferase"
)
HEADER_TAIL = ">gapit|db=tinyamr|gene=tetA|acc=|func= n/a"

EXPECTED_FASTA = (
    f"{HEADER_WRAPPED}\n{'A' * 60}\n{'A' * 10}\n"
    f"{HEADER_SHORT}\n{SEQ_SHORT}\n"
    f"{HEADER_TAIL}\n{'AC' * 30}\nG\n"
)


def make_db(tmp_path: Path, records: tuple[Record, ...] = RECORDS) -> Path:
    """A db directory holding only records.jsonl (the pre-build state)."""
    db_dir = tmp_path / "tinyamr"
    db_dir.mkdir()
    write_records(records, db_dir / "records.jsonl")
    return db_dir


def test_generate_sequences_writes_exact_deterministic_bytes(tmp_path: Path) -> None:
    """Given the 3-record fixture (metachar gene, multi-class function, empty acc),
    When generated, Then the FASTA bytes are exactly the hand-computed
    string: tagged headers, 60-column wrap, LF-only endings."""
    db_dir = make_db(tmp_path)
    generate_sequences(db_dir / "records.jsonl", db_dir / "sequences")
    assert (db_dir / "sequences").read_bytes() == EXPECTED_FASTA.encode("utf-8")


def test_generate_sequences_failure_leaves_no_partial_file(tmp_path: Path) -> None:
    """Given a valid record followed by one with a line break in its product,
    When generated, Then BUILD_INVALID names the offending gene and neither a
    partial `sequences` nor a temp file is left behind."""
    db_dir = make_db(
        tmp_path,
        (
            Record(db=DB, gene="goodGene", sequence="ACGT"),
            Record(db=DB, gene="badGene", sequence="ACGT", product="two\nlines"),
        ),
    )
    with pytest.raises(DatabaseError) as excinfo:
        generate_sequences(db_dir / "records.jsonl", db_dir / "sequences")
    error = excinfo.value
    assert error.code == "BUILD_INVALID"
    assert error.exit_code == 4
    assert error.context == {"gene": "badGene"}
    assert [entry.name for entry in db_dir.iterdir()] == ["records.jsonl"]


def test_generate_sequences_rejects_empty_records_file(tmp_path: Path) -> None:
    """Given a records file with only blank lines (zero records), When
    generated, Then BUILD_INVALID fires (makeblastdb cannot index nothing)
    with the file in context and no sequences file appears."""
    db_dir = tmp_path / "tinyamr"
    db_dir.mkdir()
    (db_dir / "records.jsonl").write_text("\n\n", encoding="utf-8")
    with pytest.raises(DatabaseError) as excinfo:
        generate_sequences(db_dir / "records.jsonl", db_dir / "sequences")
    error = excinfo.value
    assert error.code == "BUILD_INVALID"
    assert error.context == {"file": str(db_dir / "records.jsonl")}
    assert not (db_dir / "sequences").exists()


def test_generate_sequences_keeps_duplicate_db_gene_pairs(tmp_path: Path) -> None:
    """Given two records sharing (db, gene), When generated, Then both are
    kept in records.jsonl order (upstream DBs contain duplicates; no dedup)."""
    db_dir = make_db(
        tmp_path,
        (
            Record(db=DB, gene="tetA", sequence="ACGT"),
            Record(db=DB, gene="tetA", sequence="ACG"),
        ),
    )
    generate_sequences(db_dir / "records.jsonl", db_dir / "sequences")
    assert (db_dir / "sequences").read_text(encoding="utf-8") == (
        f"{HEADER_TAIL}\nACGT\n{HEADER_TAIL}\nACG\n"
    )


def test_verify_round_trip_accepts_generated_file(tmp_path: Path) -> None:
    """Given freshly generated sequences, When self-checked, Then every header
    decodes back to its record (no exception)."""
    db_dir = make_db(tmp_path)
    generate_sequences(db_dir / "records.jsonl", db_dir / "sequences")
    verify_sequences(db_dir / "sequences", db_dir / "records.jsonl", DB)


def test_verify_round_trip_flags_corrupted_escape(tmp_path: Path) -> None:
    """Given a generated file with an invalid percent escape written into one
    header, When self-checked, Then BUILD_SELF_CHECK_FAILED names the gene and
    the underlying HEADER_MALFORMED is chained."""
    db_dir = make_db(tmp_path)
    generate_sequences(db_dir / "records.jsonl", db_dir / "sequences")
    corrupted = db_dir / "sequences"
    corrupted.write_text(
        corrupted.read_text(encoding="utf-8").replace("gene=tetA", "gene=tet%ZZ"),
        encoding="utf-8",
    )
    with pytest.raises(DatabaseError) as excinfo:
        verify_sequences(corrupted, db_dir / "records.jsonl", DB)
    error = excinfo.value
    assert error.code == "BUILD_SELF_CHECK_FAILED"
    assert error.context == {"gene": "tetA", "reason": "decode:HEADER_MALFORMED"}
    assert isinstance(error.__cause__, DatabaseError)


def test_verify_round_trip_flags_rewritten_gene(tmp_path: Path) -> None:
    """Given a header whose gene value was rewritten to a DIFFERENT but valid
    gene, When self-checked, Then BUILD_SELF_CHECK_FAILED (header_mismatch),
    proving the check compares decoded fields, not just syntax."""
    db_dir = make_db(tmp_path)
    generate_sequences(db_dir / "records.jsonl", db_dir / "sequences")
    corrupted = db_dir / "sequences"
    corrupted.write_text(
        corrupted.read_text(encoding="utf-8").replace("gene=tetA", "gene=tetB"),
        encoding="utf-8",
    )
    with pytest.raises(DatabaseError) as excinfo:
        verify_sequences(corrupted, db_dir / "records.jsonl", DB)
    assert excinfo.value.context == {"gene": "tetA", "reason": "header_mismatch"}


def test_verify_round_trip_flags_mutated_sequence(tmp_path: Path) -> None:
    """Given one mutated sequence letter, When self-checked, Then
    BUILD_SELF_CHECK_FAILED (sequence_mismatch) names the record's gene."""
    db_dir = make_db(tmp_path)
    generate_sequences(db_dir / "records.jsonl", db_dir / "sequences")
    corrupted = db_dir / "sequences"
    corrupted.write_text(
        corrupted.read_text(encoding="utf-8").replace("G\n", "T\n"), encoding="utf-8"
    )
    with pytest.raises(DatabaseError) as excinfo:
        verify_sequences(corrupted, db_dir / "records.jsonl", DB)
    assert excinfo.value.context == {"gene": "tetA", "reason": "sequence_mismatch"}


def test_build_database_nucl_produces_all_artifacts(tmp_path: Path) -> None:
    """Given a 3-record db dir, When built (nucl), Then all four artifacts
    exist with the exact contract: deterministic sequences bytes, a BLAST
    index, an .mmi, and a manifest whose every field is correct and whose
    sha256 matches the sequences file."""
    db_dir = make_db(tmp_path)
    manifest = build_database(
        db_dir,
        name=DB,
        dbtype="nucl",
        source_urls=("https://example.org/tinyamr.fasta",),
        fetched_at="2026-09-17T10:30:00Z",
    )
    assert (db_dir / "sequences").read_text(encoding="utf-8") == EXPECTED_FASTA
    assert (db_dir / "sequences.nin").is_file()
    assert (db_dir / "sequences.mmi").is_file()
    on_disk = read_manifest(db_dir / "gapit-manifest.json")
    assert on_disk == manifest
    assert on_disk.schema_name == "gapit.manifest/1"
    assert on_disk.name == DB
    assert on_disk.source_urls == ("https://example.org/tinyamr.fasta",)
    assert on_disk.fetched_at == "2026-09-17T10:30:00Z"
    assert on_disk.sha256 == hashlib.sha256(EXPECTED_FASTA.encode("utf-8")).hexdigest()
    assert on_disk.n_records == 3
    assert on_disk.dbtype == "nucl"
    assert on_disk.header_format == "gapit/v1"
    assert on_disk.upstream_version == ""
    assert on_disk.tool.name == "gapit"
    assert on_disk.tool.version == __version__
    assert on_disk.makeblastdb_version
    assert on_disk.minimap2_version
    assert (
        (db_dir / "gapit-manifest.json")
        .read_text(encoding="utf-8")
        .startswith('{\n  "schema": "gapit.manifest/1",')
    )


def test_rebuild_is_byte_identical(tmp_path: Path) -> None:
    """Given the same inputs, When built twice, Then sequences and manifest
    bytes are identical (fetched_at is an input; nothing is timestamped)."""
    db_dir = make_db(tmp_path)
    first = build_database(
        db_dir, name=DB, dbtype="nucl", source_urls=(), fetched_at="2026-09-17T10:30:00Z"
    )
    sequences = (db_dir / "sequences").read_bytes()
    manifest_bytes = (db_dir / "gapit-manifest.json").read_bytes()
    second = build_database(
        db_dir, name=DB, dbtype="nucl", source_urls=(), fetched_at="2026-09-17T10:30:00Z"
    )
    assert (db_dir / "sequences").read_bytes() == sequences
    assert (db_dir / "gapit-manifest.json").read_bytes() == manifest_bytes
    assert second == first


def test_build_database_prot_skips_mmi_and_builds_pin(tmp_path: Path) -> None:
    """Given dbtype="prot" (synthetic: letters are ACGT but the manifest
    declares the type), When built, Then a .pin index exists (the mol_type
    heuristic was skipped), no .mmi is created, and the manifest records
    dbtype prot."""
    db_dir = make_db(tmp_path)
    manifest = build_database(
        db_dir, name=DB, dbtype="prot", source_urls=(), fetched_at="2026-09-17T10:30:00Z"
    )
    assert (db_dir / "sequences.pin").is_file()
    assert not (db_dir / "sequences.mmi").exists()
    assert manifest.dbtype == "prot"


def test_manifest_written_last_when_mmi_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given an .mmi build failure, When built, Then the error propagates, the
    earlier artifacts (sequences, BLAST index) exist, and NO manifest exists —
    the manifest certifies the artifacts and is written last."""
    db_dir = make_db(tmp_path)

    def fail_mmi(sequences_path: Path, mmi_path: Path) -> None:
        raise DatabaseError("simulated minimap2 crash", code="MMI_BUILD_FAILED")

    monkeypatch.setattr(dbbuild, "_build_mmi", fail_mmi)
    with pytest.raises(DatabaseError) as excinfo:
        build_database(
            db_dir, name=DB, dbtype="nucl", source_urls=(), fetched_at="2026-09-17T10:30:00Z"
        )
    assert excinfo.value.code == "MMI_BUILD_FAILED"
    assert (db_dir / "sequences").is_file()
    assert (db_dir / "sequences.nin").is_file()
    assert not (db_dir / "gapit-manifest.json").exists()
