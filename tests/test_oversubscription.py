"""Integration tests for the oversubscription note (#7): when
``--jobs * --threads`` exceeds ``os.cpu_count()``, `gapit screen` emits ONE
stderr note; it never fires when quiet or when the product stays within the
cpu count (including exactly at it, and when cpu_count is unknown)."""

import os
import shutil
from pathlib import Path

import pytest
from typer.testing import CliRunner, Result

from gapit.cli import app
from gapit.db import make_blast_db

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"
CONTIGS = Path(__file__).parent / "data" / "contigs"

runner = CliRunner()


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    """Fresh datadir with a built tinyamr index, one per test."""
    target = tmp_path / "datadir"
    shutil.copytree(FIXTURE_DB_DIR, target)
    make_blast_db(target / "tinyamr" / "sequences", "tinyamr")
    return target


def _screen_with_cpus(
    datadir: Path, monkeypatch: pytest.MonkeyPatch, cpu_count: int | None, *extra: str
) -> Result:
    """Run one screen with os.cpu_count() monkeypatched (None = unknown)."""
    monkeypatch.setattr(os, "cpu_count", lambda: cpu_count)
    return runner.invoke(app, ["screen", "--db", "tinyamr", "--datadir", str(datadir), *extra])


def test_note_fires_once_when_oversubscribed(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given 2 jobs x 4 threads on a 4-cpu machine, When screening one file,
    Then the note appears exactly once with the exact wording."""
    result = _screen_with_cpus(
        datadir,
        monkeypatch,
        4,
        "--jobs",
        "2",
        "--threads",
        "4",
        str(CONTIGS / "full.fa"),
    )
    assert result.exit_code == 0
    assert "--jobs 2 --threads 4 oversubscribes 4 cpus" in result.stderr
    assert result.stderr.count("oversubscribes") == 1


def test_note_absent_when_product_equals_cpu_count(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given 2 jobs x 4 threads on an 8-cpu machine (product == cpus), When
    screening, Then no note (the guard is strictly greater-than)."""
    result = _screen_with_cpus(
        datadir,
        monkeypatch,
        8,
        "--jobs",
        "2",
        "--threads",
        "4",
        str(CONTIGS / "full.fa"),
    )
    assert result.exit_code == 0
    assert "oversubscribes" not in result.stderr


def test_note_suppressed_by_quiet(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given an oversubscribed run with --quiet, When screening, Then the note
    is suppressed with all other stderr diagnostics."""
    result = _screen_with_cpus(
        datadir,
        monkeypatch,
        4,
        "--jobs",
        "2",
        "--threads",
        "4",
        "--quiet",
        str(CONTIGS / "full.fa"),
    )
    assert result.exit_code == 0
    assert result.stderr == ""


def test_note_absent_when_cpu_count_unknown(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given os.cpu_count() -> None (unknown), When screening an
    oversubscribed configuration, Then no note and no crash."""
    result = _screen_with_cpus(
        datadir,
        monkeypatch,
        None,
        "--jobs",
        "8",
        "--threads",
        "8",
        str(CONTIGS / "full.fa"),
    )
    assert result.exit_code == 0
    assert "oversubscribes" not in result.stderr
