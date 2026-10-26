"""Database rebuild transactions: a failed --force must not damage the old DB."""

import subprocess
import threading
from pathlib import Path

import pytest

from gapit import db as db_module
from gapit import db_build_ops
from gapit.db import discover_databases, make_blast_db
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

    with staged_database(db_dir) as build_dir:
        assert build_dir != db_dir
        assert (db_dir / "old").read_text(encoding="utf-8") == "old"
        (build_dir / "new").write_text("new", encoding="utf-8")

    assert (db_dir / "new").read_text(encoding="utf-8") == "new"
    assert not (db_dir / "old").exists()
    assert list(tmp_path.glob(".gapit-txn.backup.*")) == []
    assert list(tmp_path.glob(".gapit-txn.build.*")) == []


def test_staged_database_restores_on_failure(tmp_path: Path) -> None:
    """Given an existing database dir, When the staged write raises, Then the
    original directory is restored byte-for-byte and no backup remains."""
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    (db_dir / "old").write_text("old", encoding="utf-8")

    with pytest.raises(RuntimeError), staged_database(db_dir) as build_dir:
        (build_dir / "new").write_text("new", encoding="utf-8")
        raise RuntimeError("boom")

    assert (db_dir / "old").read_text(encoding="utf-8") == "old"
    assert not (db_dir / "new").exists()
    assert list(tmp_path.glob(".gapit-txn.backup.*")) == []


def test_staged_database_removes_a_failed_first_build(tmp_path: Path) -> None:
    """Given no existing database, When the first build raises, Then the
    partial target is removed too."""
    db_dir = tmp_path / "db"

    with pytest.raises(RuntimeError), staged_database(db_dir) as build_dir:
        (build_dir / "partial").write_text("x", encoding="utf-8")
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
    assert list(datadir.glob(".gapit-txn.backup.*")) == []


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


def test_staged_database_serializes_concurrent_writers(tmp_path: Path) -> None:
    """Given two writers for one database, When they overlap, Then the second
    waits for the first transaction to publish before staging its own build."""
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    (db_dir / "version").write_text("0", encoding="utf-8")

    first_inside = threading.Event()
    release_first = threading.Event()
    second_inside = threading.Event()
    errors: list[BaseException] = []

    def first() -> None:
        try:
            with staged_database(db_dir) as build_dir:
                first_inside.set()
                assert release_first.wait(timeout=5)
                (build_dir / "version").write_text("1", encoding="utf-8")
        except BaseException as exc:  # pragma: no cover - failure reporting
            errors.append(exc)

    def second() -> None:
        try:
            assert first_inside.wait(timeout=5)
            with staged_database(db_dir) as build_dir:
                second_inside.set()
                (build_dir / "version").write_text("2", encoding="utf-8")
        except BaseException as exc:  # pragma: no cover - failure reporting
            errors.append(exc)

    one = threading.Thread(target=first)
    two = threading.Thread(target=second)
    one.start()
    assert first_inside.wait(timeout=5)
    two.start()
    assert not second_inside.wait(timeout=0.2), "second writer entered while the lock was held"
    release_first.set()
    one.join(timeout=10)
    two.join(timeout=10)

    assert errors == []
    assert second_inside.is_set()
    assert (db_dir / "version").read_text(encoding="utf-8") == "2"


def test_discovery_ignores_hidden_transaction_directories(tmp_path: Path) -> None:
    """Given a real DB beside hidden staging/backup leftovers, When databases
    are discovered, Then only the canonical directory is listed."""
    real = tmp_path / "real"
    real.mkdir()
    (real / "sequences").write_text(">g\nACGT\n", encoding="utf-8")
    for leftover in (".gapit-txn.build.real.dead", ".gapit-txn.backup.real.dead"):
        stale = tmp_path / leftover
        stale.mkdir()
        (stale / "sequences").write_text(">stale\nACGT\n", encoding="utf-8")

    assert [database.name for database in discover_databases(tmp_path)] == ["real"]


def test_make_blast_db_keeps_the_old_index_when_makeblastdb_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given an existing index and a makeblastdb run that fails, When the
    index rebuild runs, Then the old index bytes survive (the temp prefix is
    discarded)."""
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    sequences = db_dir / "sequences"
    sequences.write_text(">s\nACGT\n", encoding="utf-8")
    old_index = db_dir / "sequences.nin"
    old_index.write_text("OLD-INDEX", encoding="utf-8")

    def fail_tool(argv: list[str]) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(argv, 1, stdout="", stderr="boom")

    monkeypatch.setattr(db_module, "run_tool", fail_tool)
    with pytest.raises(DatabaseError) as excinfo:
        make_blast_db(sequences, "db", dbtype="nucl")

    assert excinfo.value.code == "MAKEBLASTDB_FAILED"
    assert old_index.read_text(encoding="utf-8") == "OLD-INDEX"
    assert list(db_dir.glob(".gapit-index.*")) == []
