"""CLI integration tests for `gapit summary` (goldens pinned like the others)."""

import json
from datetime import UTC, datetime
from pathlib import Path

import pytest
from pydantic import TypeAdapter
from typer.testing import CliRunner, Result

from gapit.cli import app
from gapit.formats.summary import SummaryDocument, render_summary_json, render_summary_md
from gapit.summary import SummaryParams, build_summary

FIXTURES = Path(__file__).parent / "data" / "summary"
GOLDEN = Path(__file__).parent / "golden"
PINNED_NOW = datetime(2026, 9, 17, 12, 0, 0, tzinfo=UTC)
MULTI = ["sample_a.tsv", "sample_b.tsv", "empty.tsv"]

runner = CliRunner()
summary_adapter = TypeAdapter(SummaryDocument)


def summary(*extra: str) -> Result:
    return runner.invoke(app, ["summary", *extra])


@pytest.fixture(autouse=True)
def in_fixtures(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(FIXTURES)


def test_default_tsv_multi(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given several reports, When summarized with defaults, Then byte-identical
    TSV to the committed golden (itself verified against real abricate)."""
    monkeypatch.chdir(FIXTURES)
    result = summary(*MULTI)
    assert result.exit_code == 0
    assert result.stdout == (GOLDEN / "summary_multi.tsv").read_text(encoding="utf-8")


def test_dutch_tsv_single_file(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given one report, When summarized, Then dutch-mode golden TSV."""
    monkeypatch.chdir(FIXTURES)
    result = summary("multi_sample.tsv")
    assert result.exit_code == 0
    assert result.stdout == (GOLDEN / "summary_dutch.tsv").read_text(encoding="utf-8")


def test_csv_output(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --format csv on CSV reports, When summarized, Then comma-joined
    golden CSV (upstream --summary --csv shape)."""
    monkeypatch.chdir(FIXTURES)
    result = summary("--format", "csv", "sample_a.csv", "sample_b.csv")
    assert result.exit_code == 0
    assert result.stdout == (GOLDEN / "summary_multi.csv").read_text(encoding="utf-8")


def test_identity_tsv(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --identity, When summarized, Then cells are %IDENTITY values and
    the stderr note mirrors upstream."""
    monkeypatch.chdir(FIXTURES)
    result = summary("--identity", *MULTI)
    assert result.exit_code == 0
    assert "sample_a.tsv\t2\t98.75;91.00\t95.10" in result.stdout
    assert "Using %IDENTITY for the summary table instead of %COVERAGE" in result.stderr


def test_identity_note_suppressed_by_quiet() -> None:
    """Given --quiet --identity, When summarized, Then stderr is silent."""
    result = summary("--quiet", "--identity", *MULTI)
    assert result.exit_code == 0
    assert result.stderr == ""


def test_nopath_basenames_labels(tmp_path: Path) -> None:
    """Given --nopath on subdirectory inputs, When summarized, Then labels are
    basenames (row order still by full key)."""
    deep = tmp_path / "sub" / "deep_report.tsv"
    deep.parent.mkdir()
    deep.write_text((FIXTURES / "sample_a.tsv").read_text(encoding="utf-8"), encoding="utf-8")
    empty = tmp_path / "empty.tsv"
    empty.write_text((FIXTURES / "empty.tsv").read_text(encoding="utf-8"), encoding="utf-8")
    result = summary("--nopath", str(deep), str(empty))
    assert result.exit_code == 0
    labels = [line.split("\t")[0] for line in result.stdout.splitlines()[1:]]
    assert labels == ["empty.tsv", "deep_report.tsv"]


def test_duplicate_path_warns_on_stderr() -> None:
    """Given the same file twice, When summarized, Then the WARNING goes to
    stderr, the row appears once (non-dutch: keyed by input filename), and
    exit is 0 — matches the verified upstream run."""
    twice = summary("sample_a.tsv", "sample_a.tsv")
    assert twice.exit_code == 0
    assert twice.stdout == (
        "#FILE\tNUM_FOUND\tfeature_a\tfeature_b\nsample_a.tsv\t2\t99.50;52.00\t76.00\n"
    )
    assert "WARNING: Skipping duplicate file: sample_a.tsv" in twice.stderr


def test_quiet_silences_duplicate_warning() -> None:
    result = summary("--quiet", "sample_a.tsv", "sample_a.tsv")
    assert result.exit_code == 0
    assert result.stderr == ""


def test_no_files_exits_2() -> None:
    """Given no arguments, When summarized, Then USAGE_ERROR envelope, exit 2."""
    result = summary()
    assert result.exit_code == 2
    envelope = json.loads(result.stderr)
    assert envelope["code"] == "USAGE_ERROR"


def test_missing_file_exits_5() -> None:
    result = summary("nope.tsv")
    assert result.exit_code == 5
    assert json.loads(result.stderr)["code"] == "INPUT_NOT_FOUND"


def test_malformed_report_exits_5(tmp_path: Path) -> None:
    """Given a report with a short data row, When summarized, Then the
    SUMMARY_MALFORMED envelope names the file and line."""
    bad = tmp_path / "bad.tsv"
    bad.write_text(
        (FIXTURES / "empty.tsv").read_text(encoding="utf-8") + "k.fa\tone\ttwo\n",
        encoding="utf-8",
    )
    result = summary(str(bad))
    assert result.exit_code == 5
    envelope = json.loads(result.stderr)
    assert envelope["code"] == "SUMMARY_MALFORMED"
    assert envelope["context"]["file"] == str(bad)
    assert envelope["context"]["line"] == "2"


def test_json_output_validates(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --format json, When summarized, Then stdout parses as
    gapit.summary/1 with the expected matrix."""
    monkeypatch.chdir(FIXTURES)
    result = summary("--format", "json", *MULTI)
    assert result.exit_code == 0
    document = summary_adapter.validate_json(result.stdout)
    assert document.schema_name == "gapit.summary/1"
    assert document.genes == ["feature_a", "feature_b"]
    assert [row.file for row in document.rows] == ["empty.tsv", "sample_a.tsv", "sample_b.tsv"]


def test_md_output_shape(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --format md, When summarized, Then frontmatter + matrix table."""
    monkeypatch.chdir(FIXTURES)
    result = summary("--format", "md", *MULTI)
    assert result.exit_code == 0
    assert "schema: gapit.summary/1" in result.stdout
    assert "| sample_a.tsv | 2 | 99.50;52.00 | 76.00 |" in result.stdout


def test_golden_json_render(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given the pinned clock, When the CLI matrix is rendered, Then the JSON
    is byte-identical to the committed golden (CLI wiring == renderer)."""
    monkeypatch.chdir(FIXTURES)
    live = summary("--format", "json", *MULTI)
    matrix = build_summary(
        [Path(name) for name in MULTI], SummaryParams(), warn=lambda _message: None
    )
    assert live.exit_code == 0
    pinned = render_summary_json(matrix, now=PINNED_NOW)
    assert pinned == (GOLDEN / "summary_multi.json").read_text(encoding="utf-8")
    # CLI output differs from the golden only in created_at.
    assert (
        live.stdout.replace(json.loads(live.stdout)["created_at"], "2026-09-17T12:00:00Z") == pinned
    )


def test_golden_md_render(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given the pinned clock, When the matrix is rendered as MD, Then
    byte-identical to the committed golden."""
    monkeypatch.chdir(FIXTURES)
    matrix = build_summary(
        [Path(name) for name in MULTI], SummaryParams(), warn=lambda _message: None
    )
    assert render_summary_md(matrix, now=PINNED_NOW) == (GOLDEN / "summary_multi.md").read_text(
        encoding="utf-8"
    )


def test_schema_summary_registered() -> None:
    """Given `gapit schema summary`, When run, Then it prints the JSON Schema
    of gapit.summary/1."""
    result = runner.invoke(app, ["schema", "summary"])
    assert result.exit_code == 0
    schema = json.loads(result.stdout)
    assert schema["title"] == "SummaryDocument"
    assert set(schema["properties"]) >= {"schema", "tool", "created_at", "params", "genes", "rows"}


def test_schema_unknown_still_lists_summary() -> None:
    result = runner.invoke(app, ["schema", "bogus"])
    assert result.exit_code == 2
    assert "summary" in json.loads(result.stderr)["message"]
