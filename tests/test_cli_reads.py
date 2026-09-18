"""CLI tests for reads mode (gapit.reads/1) plus committed goldens."""

import json
import shutil
from datetime import UTC, datetime
from pathlib import Path

import pytest
from pydantic import TypeAdapter
from typer.testing import CliRunner

from gapit.cli import app
from gapit.formats.json import ReadsDocument, render_reads_json
from gapit.formats.md import render_reads_markdown
from gapit.reads import ReadsParams, screen_reads

READS_DB_DIR = Path(__file__).parent / "data" / "reads_db"
READS = Path(__file__).parent / "data" / "reads"
GOLDEN = Path(__file__).parent / "golden"
PINNED_NOW = datetime(2026, 9, 17, 12, 0, 0, tzinfo=UTC)

runner = CliRunner()
reads_adapter = TypeAdapter(ReadsDocument)
PARAMS = ReadsParams(db="tinyreads")


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    target = tmp_path / "datadir"
    shutil.copytree(READS_DB_DIR, target)
    return target


def envelope(stderr: str) -> dict[str, str]:
    return json.loads([line for line in stderr.splitlines() if line.strip()][-1])


def test_reads_default_json(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --r1 with no --format, When run, Then gapit.reads/1 JSON on
    stdout and Screening/Detected chatter on stderr."""
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app, ["screen", "--r1", "tetx_full.fq", "--db", "tinyreads", "--datadir", str(datadir)]
    )
    assert result.exit_code == 0
    document = reads_adapter.validate_json(result.stdout)
    assert document.schema_name == "gapit.reads/1"
    assert document.files[0].reads == ["tetx_full.fq"]
    (entry,) = document.files[0].genes
    assert entry.gene == "tetX"
    assert entry.present is True
    assert "Screening reads:" in result.stderr
    assert "Detected 1 present genes" in result.stderr


def test_reads_and_positional_files_are_mutually_exclusive(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "tetx_full.fq",
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
            "junk.fq",
        ],
    )
    assert result.exit_code == 2
    assert envelope(result.stderr)["code"] == "USAGE_ERROR"


def test_reads_mode_rejects_tsv_and_csv(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(READS)
    for fmt in ("tsv", "csv"):
        result = runner.invoke(
            app,
            [
                "screen",
                "--r1",
                "tetx_full.fq",
                "--db",
                "tinyreads",
                "--datadir",
                str(datadir),
                "--format",
                fmt,
            ],
        )
        assert result.exit_code == 2
        assert envelope(result.stderr)["code"] == "USAGE_ERROR"


def test_missing_reads_file_exits_5(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        ["screen", "--r1", "nope.fq", "--db", "tinyreads", "--datadir", str(datadir)],
    )
    assert result.exit_code == 5
    assert envelope(result.stderr)["code"] == "INPUT_NOT_FOUND"


def test_reads_md_output(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "tetx_full.fq",
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
            "--format",
            "md",
        ],
    )
    assert result.exit_code == 0
    assert result.stdout.startswith("---\n")
    assert "| Gene | Breadth% | Depth | Reads | Present |" in result.stdout


def test_r2_without_r1_exits_2(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        ["screen", "--r2", "tetx_R2.fq", "--db", "tinyreads", "--datadir", str(datadir)],
    )
    assert result.exit_code == 2
    assert envelope(result.stderr)["code"] == "USAGE_ERROR"


def test_r2_lane_count_mismatch_exits_2(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "tetx_lane1.fq,tetx_lane2.fq",
            "--r2",
            "tetx_R2.fq",
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
        ],
    )
    assert result.exit_code == 2
    assert envelope(result.stderr)["code"] == "USAGE_ERROR"


def test_empty_comma_element_exits_2(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "tetx_lane1.fq,,tetx_lane2.fq",
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
        ],
    )
    assert result.exit_code == 2
    assert envelope(result.stderr)["code"] == "USAGE_ERROR"


def test_multi_lane_se_union(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given two comma-listed single-end lanes that are each sub-threshold,
    When screened as one sample, Then one gene entry whose union crosses the
    threshold and whose reads_mapped counts both lanes."""
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "tetx_lane1.fq,tetx_lane2.fq",
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
        ],
    )
    assert result.exit_code == 0
    document = reads_adapter.validate_json(result.stdout)
    (entry,) = document.files[0].genes
    assert entry.gene == "tetX"
    assert entry.present is True
    assert entry.reads_mapped == 12
    assert document.files[0].reads == ["tetx_lane1.fq", "tetx_lane2.fq"]


def test_multi_lane_pe_union(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given two paired lanes (2x2 files), When screened, Then one gene entry
    found once with all four files' reads counted."""
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "tetx_pe1_R1.fq,tetx_pe2_R1.fq",
            "--r2",
            "tetx_pe1_R2.fq,tetx_pe2_R2.fq",
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
        ],
    )
    assert result.exit_code == 0
    document = reads_adapter.validate_json(result.stdout)
    (entry,) = document.files[0].genes
    assert entry.present is True
    assert entry.reads_mapped == 12
    assert document.files[0].reads == [
        "tetx_pe1_R1.fq",
        "tetx_pe2_R1.fq",
        "tetx_pe1_R2.fq",
        "tetx_pe2_R2.fq",
    ]


def test_reads_flag_is_gone_from_help() -> None:
    """Given --help for screen, When inspected, Then --reads is absent and
    --r1/--r2 are present."""
    result = runner.invoke(app, ["screen", "--help"])
    assert result.exit_code == 0
    assert "--reads" not in result.stdout
    assert "--r1" in result.stdout
    assert "--r2" in result.stdout


def test_reads_debug_echoes_minimap2_argv(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --r1 with --debug, When run, Then the minimap2 argv is echoed
    to stderr as a `gapit: run:` line and stdout stays identical to the
    plain run (debug is stderr-only); the plain run emits no argv lines."""
    args = ["screen", "--r1", "tetx_full.fq", "--db", "tinyreads", "--datadir", str(datadir)]
    monkeypatch.chdir(READS)
    debug_run = runner.invoke(app, [*args, "--debug"])
    plain_run = runner.invoke(app, args)
    assert debug_run.exit_code == 0
    run_lines = [line for line in debug_run.stderr.splitlines() if line.startswith("gapit: run:")]
    assert run_lines and run_lines[0].startswith("gapit: run: minimap2 -x sr ")
    assert debug_run.stdout == plain_run.stdout
    assert "gapit: run:" not in plain_run.stderr


def test_golden_reads_json(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Given the tetx fixture run and a pinned now, When rendered, Then the
    JSON is byte-identical to the committed golden."""
    target = tmp_path / "datadir"
    shutil.copytree(READS_DB_DIR, target)
    from gapit.db import discover_databases

    (database,) = discover_databases(target)
    monkeypatch.chdir(READS)
    report = screen_reads(
        [(Path("tetx_full.fq"), None)], database, read_type="sr", min_breadth=90.0, threads=1
    )
    output = render_reads_json([report], PARAMS, now=PINNED_NOW)
    assert output == (GOLDEN / "reads_tinyamr.json").read_text(encoding="utf-8")
    assert reads_adapter.validate_json(output).schema_name == "gapit.reads/1"


def test_golden_reads_md(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    target = tmp_path / "datadir"
    shutil.copytree(READS_DB_DIR, target)
    from gapit.db import discover_databases

    (database,) = discover_databases(target)
    monkeypatch.chdir(READS)
    report = screen_reads(
        [(Path("tetx_full.fq"), None)], database, read_type="sr", min_breadth=90.0, threads=1
    )
    output = render_reads_markdown([report], PARAMS, now=PINNED_NOW)
    assert output == (GOLDEN / "reads_tinyamr.md").read_text(encoding="utf-8")
