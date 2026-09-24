"""Unit tests for the reads use-case perform contract: run_screen_reads
returns the rendered document instead of echoing it (MCP groundwork)."""

import shutil
from datetime import UTC, datetime
from pathlib import Path

import pytest

from gapit import screening_reads
from gapit.screening_reads import run_screen_reads

READS_DB_DIR = Path(__file__).parent / "data" / "reads_db"
READS = Path(__file__).parent / "data" / "reads"
GOLDEN = Path(__file__).parent / "golden"
PINNED_NOW = datetime(2026, 9, 17, 12, 0, 0, tzinfo=UTC)


class _FrozenDatetime:
    """Stand-in for the datetime class as _screen_lanes uses it: now(tz)."""

    @staticmethod
    def now(tz: object) -> datetime:
        return PINNED_NOW


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    target = tmp_path / "datadir"
    shutil.copytree(READS_DB_DIR, target)
    return target


def test_run_screen_reads_returns_golden_json(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the tetx fixture and a pinned clock, When run_screen_reads is
    called directly (no CliRunner) with already-split lane lists, Then it
    returns the golden JSON string for the caller to echo."""
    monkeypatch.setattr(screening_reads, "datetime", _FrozenDatetime)
    monkeypatch.chdir(READS)
    # quiet=True only silences stderr diagnostics; the returned document is
    # byte-identical either way, and direct calls must not spam real stderr.
    output = run_screen_reads(
        [Path("tetx_full.fq")],
        None,
        "tinyreads",
        datadir,
        None,
        90.0,
        0.0,
        0,
        1,
        None,
        True,
    )
    assert output == (GOLDEN / "reads_tinyamr.json").read_text(encoding="utf-8")
