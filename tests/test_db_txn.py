"""Database rebuild transactions: a failed --force must not damage the old DB."""

from pathlib import Path

import pytest

from gapit import db as db_module
from gapit import db_build_ops
from gapit.db import make_blast_db
from gapit.db_build_ops import perform_build
from gapit.db_txn import staged_database
from gapit.errors import DatabaseError, DependencyError

SEQ = "ACGT" * 20


def test_staged_database_replaces_on_success(tmp_path: Path) -> None:
    """Given an existing database dir, When a staged write succeeds, Then the
    new contents win and no backup directory remains."""
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    (db_dir / "old").write_text("old", encoding="utf-8")

    with staged_database(db_dir):
        (db_dir / "new").write_text("new", encoding="utf-8")

    assert (db_dir / "new").read_text(encoding="utf-8") == "new"
    assert not (db_dir / "old").exists()
    assert list(tmp_path.glob(".db.backup.*")) == []


def test_staged_database_restores_on_failure(tmp_path: Path) -> None:
    """Given an existing database dir, When the staged write raises, Then the
    original directory is restored byte-for-byte and no backup remains."""
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    (db_dir / "old").write_text("old", encoding="utf-8")

    with pytest.raises(RuntimeError), staged_database(db_dir):
        (db_dir / "new").write_text("new", encoding="utf-8")
        raise RuntimeError("boom")

    assert (db_dir / "old").read_text(encoding="utf-8") == "old"
    assert not (db_dir / "new").exists()
    assert list(tmp_path.glob(".db.backup.*")) == []


def test_staged_database_removes_a_failed_first_build(tmp_path: Path) -> None:
    """Given no existing database, When the first build raises, Then the
    partial target is removed too."""
    db_dir = tmp_path / "db"

    with pytest.raises(RuntimeError), staged_database(db_dir):
        (db_dir / "partial").write_text("x", encoding="utf-8")
        raise RuntimeError("boom")

    assert not db_dir.exists()


def test_force_rebuild_failure_restores_the_previous_database(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a built database, When --force rebuild fails mid-pipeline, Then
    the old sequences/manifest survive and the receipt never claims success."""
    datadir = tmp_path / "dd"
    old = datadir / "keep"
    old.mkdir(parents=True)
    (old / "sequences").write_text("OLD-SEQUENCES", encoding="utf-8")
    (old / "gapit-manifest.json").write_text("OLD-MANIFEST", encoding="utf-8")
    fasta = tmp_path / "new.fa"
    fasta.write_text(f">gapit|db=x|gene=g1|acc=A|func=f\n{SEQ}\n", encoding="utf-8")

    def fail_build(*args: object, **kwargs: object) -> object:
        raise DatabaseError("simulated build failure", code="BUILD_SELF_CHECK_FAILED")

    monkeypatch.setattr(db_build_ops, "build_database", fail_build)
    with pytest.raises(DatabaseError):
        perform_build(
            "keep",
            fasta,
            None,
            None,
            "",
            datadir,
            True,
            warn=lambda message: None,
            quiet=True,
        )

    assert (old / "sequences").read_text(encoding="utf-8") == "OLD-SEQUENCES"
    assert (old / "gapit-manifest.json").read_text(encoding="utf-8") == "OLD-MANIFEST"
    assert list(datadir.glob(".keep.backup.*")) == []


def test_make_blast_db_does_not_delete_an_index_before_preflight(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given no makeblastdb on PATH, When an index rebuild runs, Then the
    preflight fails before the existing index files are unlinked."""
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    sequences = db_dir / "sequences"
    sequences.write_text(">s\nACGT\n", encoding="utf-8")
    index = db_dir / "sequences.nin"
    index.write_text("OLD-INDEX", encoding="utf-8")

    def no_tool(name: str) -> str | None:
        return None

    def fail_tool(argv: list[str]) -> object:
        raise DependencyError(
            f"required binary not found on PATH: {argv[0]}",
            code="MISSING_DEPENDENCY",
        )

    monkeypatch.setattr(db_module.shutil, "which", no_tool)
    monkeypatch.setattr(db_module, "run_tool", fail_tool)
    with pytest.raises(DependencyError):
        make_blast_db(sequences, "db", dbtype="nucl")

    assert index.read_text(encoding="utf-8") == "OLD-INDEX"
