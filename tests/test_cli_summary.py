"""CLI integration tests for `gapit summary` (goldens pinned like the others)."""

import json
from datetime import UTC, datetime
from pathlib import Path

import pytest
from pydantic import TypeAdapter
from typer.testing import CliRunner, Result

from gapit import cmd_summary
from gapit.cli import app
from gapit.formats.summary import SummaryDocument, render_summary_json, render_summary_md
from gapit.summary import SummaryParams, build_summary

FIXTURES = Path(__file__).parent / "data" / "summary"
TYPING_DATA = Path(__file__).parent / "data" / "typing"
GOLDEN = Path(__file__).parent / "golden"
PINNED_NOW = datetime(2026, 9, 17, 12, 0, 0, tzinfo=UTC)
MULTI = ["sample_a.tsv", "sample_b.tsv", "empty.tsv"]

runner = CliRunner()
summary_adapter = TypeAdapter(SummaryDocument)


def summary(*extra: str, input: str | None = None) -> Result:
    return runner.invoke(app, ["summary", *extra], input=input)


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


def test_combined_identity_coverage_cells(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given -ic (both cell metrics), When summarized, Then each hit renders
    "identity/coverage" pairs ';'-joined, with the stderr note."""
    monkeypatch.chdir(FIXTURES)
    result = summary("-ic", *MULTI)
    assert result.exit_code == 0
    assert "Using %IDENTITY/%COVERAGE per hit" in result.stderr
    assert "sample_a.tsv\t2\t98.75/99.50;91.00/52.00\t95.10/76.00" in result.stdout
    assert "sample_b.tsv\t2\t97.00/90.00\t99.99/100.00" in result.stdout


def test_combined_short_flag_bundle_equals_both_flags(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given the Unix bundle -ic, When compared with -i -c, Then outputs
    are byte-identical."""
    monkeypatch.chdir(FIXTURES)
    bundled = summary("-ic", *MULTI)
    split = summary("-i", "-c", *MULTI)
    assert bundled.exit_code == 0 and split.exit_code == 0
    assert bundled.stdout == split.stdout


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
    assert "Using %IDENTITY for the summary table instead of +/- presence" in result.stderr


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
    assert twice.stdout == ("#FILE\tNUM_FOUND\tfeature_a\tfeature_b\nsample_a.tsv\t2\t+\t+\n")
    assert "WARNING: Skipping duplicate file: sample_a.tsv" in twice.stderr


def test_quiet_silences_duplicate_warning() -> None:
    result = summary("--quiet", "sample_a.tsv", "sample_a.tsv")
    assert result.exit_code == 0
    assert result.stderr == ""


def test_no_files_with_flag_reads_stdin() -> None:
    """Given a flag but no report files (an empty pipe on stdin), When
    summarized, Then the empty matrix renders — an omitted input now means
    stdin, not a usage error (the `gapit typing` precedent)."""
    result = summary("--quiet", input="")
    assert result.exit_code == 0
    assert result.stdout == "#FILE\tNUM_FOUND\n"


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
    assert "| sample_a.tsv | 2 | + | + |" in result.stdout


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


class TestStdinInput:
    """The piped report table: `gapit screen -d ecoli_dec *.fna | gapit
    summary` summarizes stdin with no file arguments; `-` is the explicit
    marker (it reads stdin even under a terminal, until EOF)."""

    def test_bare_invocation_at_a_terminal_prints_help(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Given no arguments at all with a terminal on stdin, When summary
        runs, Then the full help prints and exits 2 with no error envelope
        (the help-on-bare policy survives the pipe support; the tty probe is
        pinned because a CliRunner stdin is never a tty)."""
        monkeypatch.setattr(cmd_summary, "_stdin_is_tty", lambda: True)
        result = runner.invoke(app, ["summary"])
        assert result.exit_code == 2
        assert "Usage" in result.stdout
        assert "gapit.error" not in result.stderr

    def test_piped_batched_table_summarizes_per_file(self) -> None:
        """Given a batched multi-FILE screen table on stdin (one input),
        When summarized, Then dutch mode applies — one row per FILE value,
        byte-identical to summarizing the same table by path."""
        table = (FIXTURES / "multi_sample.tsv").read_text(encoding="utf-8")
        via_file = summary("multi_sample.tsv")
        piped = summary(input=table)
        assert piped.exit_code == 0, piped.stderr
        assert piped.stdout == via_file.stdout
        assert piped.stdout == (GOLDEN / "summary_dutch.tsv").read_text(encoding="utf-8")

    def test_explicit_dash_matches_the_file_run(self) -> None:
        """Given the `-` marker, When summary runs, Then stdin is the table
        and the output matches the file-based run byte for byte."""
        table = (FIXTURES / "multi_sample.tsv").read_text(encoding="utf-8")
        expected = summary("multi_sample.tsv")
        dash = summary("-", input=table)
        assert dash.exit_code == 0, dash.stderr
        assert dash.stdout == expected.stdout

    def test_dash_with_a_terminal_still_reads(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Given `-` while stdin is a terminal, When summary runs, Then the
        marker still reads stdin until EOF — the tty guard only covers the
        argument-less invocation."""
        monkeypatch.setattr(cmd_summary, "_stdin_is_tty", lambda: True)
        table = (FIXTURES / "multi_sample.tsv").read_text(encoding="utf-8")
        result = summary("-", input=table)
        assert result.exit_code == 0, result.stderr
        assert result.stdout.startswith("#FILE\tNUM_FOUND")

    def test_mixed_dash_and_files_is_a_usage_error(self) -> None:
        """Given `-` alongside a file argument, When summary runs, Then a
        usage error demands one or the other (v1: either stdin or files,
        never both)."""
        result = summary("-", "sample_a.tsv", input="")
        assert result.exit_code == 2
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "USAGE_ERROR"
        assert "not both" in envelope["message"]

    def test_empty_stdin_is_the_empty_matrix(self) -> None:
        """Given an empty pipe, When summarized, Then the empty matrix
        renders (header only) — the same output an empty table file gives."""
        result = summary(input="")
        assert result.exit_code == 0
        assert result.stdout == "#FILE\tNUM_FOUND\n"

    def test_malformed_piped_row_names_stdin(self) -> None:
        """Given a malformed piped row, When summarized, Then the
        SUMMARY_MALFORMED envelope names `-` and the line number."""
        table = (FIXTURES / "empty.tsv").read_text(encoding="utf-8") + "k.fa\tone\ttwo\n"
        result = summary(input=table)
        assert result.exit_code == 5
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "SUMMARY_MALFORMED"
        assert "(-:2)" in envelope["message"]
        assert envelope["context"]["file"] == "-"

    def test_screen_pipe_summarizes_identically_to_the_table(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Given one real screen run over two isolates delivered three
        ways — written with `-o` then summarized by path, summarized from
        the redirected table on stdin, and piped bare (`screen | summary`)
        — When compared, Then all three matrices are byte-identical TSV
        (the canonical rightsholder pipe)."""
        monkeypatch.chdir(tmp_path)
        datadir = tmp_path / "datadir"
        datadir.mkdir()
        samples = [
            str(TYPING_DATA / "dec_s3_stx2a_escV_aggR_uidA.fasta"),
            str(TYPING_DATA / "dec_s2_pic_astA_uidA.fasta"),
        ]
        table_path = tmp_path / "combined.tsv"
        to_file = runner.invoke(
            app,
            [
                "screen",
                *samples,
                "--db",
                "ecoli_dec",
                "--datadir",
                str(datadir),
                "--nopath",
                "-o",
                str(table_path),
                "--quiet",
            ],
        )
        assert to_file.exit_code == 0, to_file.stderr
        on_stdout = runner.invoke(
            app,
            ["screen", *samples, "--db", "ecoli_dec", "--datadir", str(datadir), "--nopath", "-q"],
        )
        assert on_stdout.exit_code == 0, on_stdout.stderr

        by_path = summary(str(table_path))
        redirected = summary(input=table_path.read_text(encoding="utf-8"))
        piped = summary(input=on_stdout.stdout)
        assert by_path.exit_code == 0, by_path.stderr
        assert redirected.exit_code == 0, redirected.stderr
        assert piped.exit_code == 0, piped.stderr
        assert by_path.stdout == redirected.stdout == piped.stdout
        assert by_path.stdout.startswith("#FILE\tNUM_FOUND\t")
        assert {line.split("\t")[0] for line in by_path.stdout.splitlines()[1:]} == {
            "dec_s2_pic_astA_uidA.fasta",
            "dec_s3_stx2a_escV_aggR_uidA.fasta",
        }
