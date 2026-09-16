"""Tests for the Markdown screening report."""

import shutil
from datetime import UTC, datetime
from pathlib import Path

import pytest

from gapit.blast import screen_file
from gapit.db import Database, make_blast_db
from gapit.formats.md import render_markdown
from gapit.report import Report, ScreeningParams

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"
CONTIGS = Path(__file__).parent / "data" / "contigs"
GOLDEN = Path(__file__).parent / "golden"
MULTI_FILES = ("full.fa", "gap.fa", "none.fa", "sort.fa")
PINNED_NOW = datetime(2026, 9, 17, 12, 0, 0, tzinfo=UTC)
PARAMS = ScreeningParams(db="tinyamr")

HEADER_ROW = (
    "| Sequence | Start | End | Strand | Gene | Coverage | Map | Gaps | %Coverage | "
    "%Identity | Database | Accession | Product | Resistance |"
)


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    target = tmp_path / "datadir"
    shutil.copytree(FIXTURE_DB_DIR, target)
    make_blast_db(target / "tinyamr" / "sequences", "tinyamr")
    return target


@pytest.fixture()
def reports(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> list[Report]:
    monkeypatch.chdir(CONTIGS)
    database = Database(
        name="tinyamr", path=datadir / "tinyamr", sequences_path=datadir / "tinyamr" / "sequences"
    )
    return [screen_file(Path(name), database, PARAMS) for name in MULTI_FILES]


def frontmatter(output: str) -> dict[str, str]:
    lines = output.splitlines()
    assert lines[0] == "---"
    end = lines.index("---", 1)
    entries: dict[str, str] = {}
    for line in lines[1:end]:
        key, _, value = line.partition(": ")
        entries[key] = value
    return entries


def test_golden_markdown(reports: list[Report]) -> None:
    """Given the multi-file fixture run and a pinned now, When rendered, Then
    the output is byte-identical to the committed golden Markdown."""
    assert render_markdown(reports, PARAMS, now=PINNED_NOW) == (
        GOLDEN / "screen_multi.md"
    ).read_text(encoding="utf-8")


def test_frontmatter_keys_and_values(reports: list[Report]) -> None:
    output = render_markdown(reports, PARAMS, now=PINNED_NOW)
    meta = frontmatter(output)
    assert list(meta) == [
        "schema",
        "tool",
        "created_at",
        "db",
        "minid",
        "mincov",
        "threads",
        "files",
        "hits",
    ]
    assert meta["schema"] == "gapit.report/1"
    assert meta["tool"].startswith("gapit ")
    assert meta["created_at"] == "2026-09-17T12:00:00Z"
    assert meta["db"] == "tinyamr"
    assert meta["files"] == "4"
    assert meta["hits"] == "6"


def test_table_header_row_is_exact(reports: list[Report]) -> None:
    output = render_markdown(reports, PARAMS, now=PINNED_NOW)
    assert HEADER_ROW in output.splitlines()
    assert "# gapit screening report" in output.splitlines()
    assert "## `sort.fa`" in output.splitlines()


def test_zero_hit_file_says_no_hits(reports: list[Report]) -> None:
    output = render_markdown(reports, PARAMS, now=PINNED_NOW)
    lines = output.splitlines()
    assert "_No hits._" in lines
    index = lines.index("## `none.fa`")
    assert lines[index + 2] == "_No hits._"


def test_lf_endings_and_trailing_newline(reports: list[Report]) -> None:
    output = render_markdown(reports, PARAMS, now=PINNED_NOW)
    assert "\r" not in output
    assert output.endswith("\n")


def test_cli_format_md(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from typer.testing import CliRunner

    from gapit.cli import app

    monkeypatch.chdir(CONTIGS)
    result = CliRunner().invoke(
        app,
        ["screen", *MULTI_FILES, "--db", "tinyamr", "--datadir", str(datadir), "--format", "md"],
    )
    assert result.exit_code == 0
    assert result.stdout.startswith("---\n")
    assert HEADER_ROW in result.stdout
