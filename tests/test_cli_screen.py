"""Integration tests: gapit screen CLI over the real pipeline on tinyamr fixtures."""

import gzip
import json
import re
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import Literal

import pytest
from pydantic import TypeAdapter
from typer.testing import CliRunner, Result

from gapit import screening as gapit_screening
from gapit.cli import app
from gapit.db import Database, make_blast_db
from gapit.errors import ErrorEnvelope, InputError
from gapit.report import Report, ScreeningParams
from gapit.screening import OutputFormat, run_screen

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"
CONTIGS = Path(__file__).parent / "data" / "contigs"
GOLDEN = Path(__file__).parent / "golden"
MULTI_FILES = ["full.fa", "gap.fa", "none.fa", "sort.fa"]

runner = CliRunner()
envelope_adapter = TypeAdapter(ErrorEnvelope)
ANSI_STYLE = re.compile(r"\x1b\[[0-9;]*m")


def last_envelope(stderr: str) -> ErrorEnvelope:
    """Parse the last non-empty stderr line as the gapit.error/1 envelope."""
    lines = [line for line in stderr.splitlines() if line.strip()]
    assert lines, "expected an error envelope on stderr"
    return envelope_adapter.validate_json(lines[-1])


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    """Fresh datadir with a built tinyamr index, one per test."""
    target = tmp_path / "datadir"
    shutil.copytree(FIXTURE_DB_DIR, target)
    make_blast_db(target / "tinyamr" / "sequences", "tinyamr")
    return target


def screen(datadir: Path, *extra: str) -> Result:
    """Invoke `gapit screen ... --db tinyamr --datadir <dd> [extra]`."""
    return runner.invoke(app, ["screen", "--db", "tinyamr", "--datadir", str(datadir), *extra])


def test_default_tsv_run(datadir: Path) -> None:
    """Given sort.fa, When screened with defaults, Then exit 0, header + 4 rows,
    Processing/Found chatter on stderr."""
    result = screen(datadir, str(CONTIGS / "sort.fa"))
    assert result.exit_code == 0
    lines = result.stdout.splitlines()
    assert lines[0].startswith("#FILE\tSEQUENCE\t")
    assert len(lines) == 5
    assert lines[1].startswith(f"{CONTIGS / 'sort.fa'}\tcontigA\t1\t79\t+\ttetA\t")
    assert "Processing:" in result.stderr
    assert "Found 4 genes" in result.stderr


def test_golden_tsv_multi_file_nopath(datadir: Path) -> None:
    """Given all four contig fixtures with --nopath, When screened, Then stdout
    is byte-identical to the committed golden TSV."""
    result = screen(datadir, "--nopath", *[str(CONTIGS / name) for name in MULTI_FILES])
    assert result.exit_code == 0
    assert result.stdout == (GOLDEN / "tinyamr_multi_nopath.tsv").read_text(encoding="utf-8")


def test_golden_csv_multi_file_nopath(datadir: Path) -> None:
    """Given the same run with --format csv, When screened, Then stdout is
    byte-identical to the committed golden CSV."""
    result = screen(
        datadir, "--format", "csv", "--nopath", *[str(CONTIGS / name) for name in MULTI_FILES]
    )
    assert result.exit_code == 0
    assert result.stdout == (GOLDEN / "tinyamr_multi_nopath.csv").read_text(encoding="utf-8")


def test_run_screen_returns_cli_identical_bytes(datadir: Path) -> None:
    """Given a direct run_screen call (the MCP route), When compared with the
    CLI twin run, Then the returned string is byte-identical to CLI stdout —
    the perform contract: the use-case renders, the caller echoes."""
    result = screen(datadir, "--nopath", *[str(CONTIGS / name) for name in MULTI_FILES])
    assert result.exit_code == 0
    output = run_screen(
        [CONTIGS / name for name in MULTI_FILES],
        "tinyamr",
        datadir,
        80.0,
        80.0,
        1,
        1,
        None,
        quiet=True,
        noheader=False,
        nopath=True,
        debug=False,
        output_format=OutputFormat.tsv,
    )
    assert output == result.stdout


def test_csv_flag_is_gone(datadir: Path) -> None:
    """Given the removed --csv flag, When screened, Then it is rejected as an
    unknown option with exit 2 and no data on stdout (lock-in: --format csv is
    the only spelling)."""
    result = screen(datadir, "--csv", "--nopath", str(CONTIGS / "full.fa"))
    assert result.exit_code == 2
    assert result.stdout == ""


def test_noheader_suppresses_header_line(datadir: Path) -> None:
    result = screen(datadir, "--noheader", str(CONTIGS / "full.fa"))
    assert result.exit_code == 0
    assert "#FILE" not in result.stdout
    assert len(result.stdout.splitlines()) == 1


def test_nopath_basenames_file_column(datadir: Path) -> None:
    result = screen(datadir, "--nopath", str(CONTIGS / "full.fa"))
    assert result.stdout.splitlines()[1].startswith("full.fa\t")


def test_fofn_replaces_positional_files(datadir: Path, tmp_path: Path) -> None:
    """Given a --fofn listing full.fa AND a positional sort.fa, When screened,
    Then only the fofn files run (upstream behavior)."""
    fofn = tmp_path / "files.txt"
    fofn.write_text(f"{CONTIGS / 'full.fa'}\n\n{CONTIGS / 'none.fa'}\n", encoding="utf-8")
    result = screen(datadir, "--nopath", "--fofn", str(fofn), str(CONTIGS / "sort.fa"))
    assert result.exit_code == 0
    files_in_output = {line.split("\t")[0] for line in result.stdout.splitlines()[1:]}
    assert files_in_output == {"full.fa"}
    assert result.stdout.count("#FILE") == 1


def test_quiet_silences_stderr_diagnostics(datadir: Path) -> None:
    result = screen(datadir, "--quiet", str(CONTIGS / "full.fa"))
    assert result.exit_code == 0
    assert result.stderr == ""


def test_invalid_minid_exits_2(datadir: Path) -> None:
    result = screen(datadir, "--minid", "0", str(CONTIGS / "full.fa"))
    assert result.exit_code == 2
    assert json.loads(result.stderr)["code"] == "USAGE_ERROR"


def test_invalid_mincov_exits_2(datadir: Path) -> None:
    result = screen(datadir, "--mincov", "101", str(CONTIGS / "full.fa"))
    assert result.exit_code == 2


def test_invalid_threads_exits_2(datadir: Path) -> None:
    result = screen(datadir, "--threads", "0", str(CONTIGS / "full.fa"))
    assert result.exit_code == 2


def test_missing_input_file_exits_5(datadir: Path) -> None:
    result = screen(datadir, str(datadir / "nope.fa"))
    assert result.exit_code == 5
    assert json.loads(result.stderr)["code"] == "INPUT_NOT_FOUND"


def test_junk_input_exits_5(datadir: Path, tmp_path: Path) -> None:
    junk = tmp_path / "junk.txt"
    junk.write_text("this is not sequence data at all\n", encoding="utf-8")
    result = screen(datadir, str(junk))
    assert result.exit_code == 5


def test_truncated_gz_contigs_exits_5_invalid_input(datadir: Path, tmp_path: Path) -> None:
    """Given a contigs .gz cut to ~60% of its bytes (partial download), When
    screened, Then exit 5 with envelope code INVALID_INPUT — the typed input
    contract, not an uncaught mid-read EOFError surfacing as UNEXPECTED/1."""
    blob = gzip.compress((CONTIGS / "full.fa").read_bytes())
    truncated = tmp_path / "truncated.fa.gz"
    truncated.write_bytes(blob[: len(blob) * 3 // 5])
    result = screen(datadir, str(truncated))
    assert result.exit_code == 5
    assert last_envelope(result.stderr).code == "INVALID_INPUT"


def test_unknown_db_exits_4_and_lists_available(datadir: Path) -> None:
    result = screen(datadir, "--db", "nope", str(CONTIGS / "full.fa"))
    assert result.exit_code == 4
    envelope = json.loads(result.stderr)
    assert envelope["code"] == "DATABASE_NOT_FOUND"
    assert "tinyamr" in envelope["message"]
    assert "Available" in envelope["message"]


def test_no_input_files_exits_2(datadir: Path) -> None:
    result = screen(datadir)
    assert result.exit_code == 2
    assert json.loads(result.stderr)["code"] == "USAGE_ERROR"


def test_missing_db_option_is_usage_error_exit_2(tmp_path: Path) -> None:
    """Given a screen invocation without --db (no default since the
    breaking change), When run, Then typer's missing-option usage error
    exits 2 and no data reaches stdout."""
    result = runner.invoke(app, ["screen", str(CONTIGS / "full.fa"), "--datadir", str(tmp_path)])
    assert result.exit_code == 2
    combined = result.stdout + result.stderr
    assert "--db" in combined
    assert result.stdout == ""


def test_help_lists_output_flag_and_required_db() -> None:
    """Given `gapit screen --help`, When inspected, Then --output is
    documented and --db is marked required (click renders [required])."""
    result = runner.invoke(app, ["screen", "--help"], env={"COLUMNS": "100"})
    assert result.exit_code == 0
    stripped = ANSI_STYLE.sub("", result.stdout)
    assert "--output" in stripped
    assert "--db" in stripped
    assert "[required]" in stripped


def test_output_file_receives_data_stdout_stays_empty(datadir: Path, tmp_path: Path) -> None:
    """Given --output PATH on a multi-file tsv run, When screened, Then the
    file holds exactly the bytes a plain run prints and stdout holds NO
    data (stdout purity with --output; stderr chatter unchanged)."""
    out_path = tmp_path / "out.tsv"
    files = [str(CONTIGS / name) for name in MULTI_FILES]
    plain = screen(datadir, "--nopath", *files)
    assert plain.exit_code == 0
    result = runner.invoke(
        app,
        [
            "screen",
            "--db",
            "tinyamr",
            "--datadir",
            str(datadir),
            "--nopath",
            "--output",
            str(out_path),
            *files,
        ],
    )
    assert result.exit_code == 0
    assert result.stdout == ""
    assert "Processing:" in result.stderr
    assert out_path.read_text(encoding="utf-8") == plain.stdout


def test_output_json_writes_single_document(datadir: Path, tmp_path: Path) -> None:
    """Given --output with --format json, When screened, Then the file is
    the single buffered gapit.report/1 document and stdout is empty."""
    out_path = tmp_path / "out.json"
    result = runner.invoke(
        app,
        [
            "screen",
            str(CONTIGS / "full.fa"),
            "--db",
            "tinyamr",
            "--datadir",
            str(datadir),
            "--quiet",
            "--format",
            "json",
            "--output",
            str(out_path),
        ],
    )
    assert result.exit_code == 0
    assert result.stdout == ""
    document = json.loads(out_path.read_text(encoding="utf-8"))
    assert document["schema"] == "gapit.report/1"
    assert document["params"]["db"] == "tinyamr"


def test_output_md_persists_streamed_prefix_on_midbatch_failure(
    datadir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given --output --format md where file 2 fails mid-batch, When
    screened, Then the file keeps the static frontmatter plus file 1's
    already-streamed section while stdout stays empty — md writes are
    incremental per chunk, exactly like the tsv path."""
    out_path = tmp_path / "out.md"
    real = gapit_screening.screen_file

    def fail_on_second(
        query: Path,
        database: Database,
        params: ScreeningParams,
        *,
        dbtype: Literal["nucl", "prot"],
        debug: bool = False,
        merge_fragments: bool = False,
    ) -> Report:
        if query.name == "gap.fa":
            raise InputError(f"boom in {query}", code="INVALID_INPUT", context={})
        return real(
            query, database, params, dbtype=dbtype, debug=debug, merge_fragments=merge_fragments
        )

    monkeypatch.setattr(gapit_screening, "screen_file", fail_on_second)
    result = runner.invoke(
        app,
        [
            "screen",
            str(CONTIGS / "full.fa"),
            str(CONTIGS / "gap.fa"),
            "--db",
            "tinyamr",
            "--datadir",
            str(datadir),
            "--nopath",
            "--format",
            "md",
            "--output",
            str(out_path),
        ],
    )
    assert result.exit_code == 5
    assert result.stdout == ""
    partial = out_path.read_text(encoding="utf-8")
    assert partial.startswith("---\nschema: gapit.report/1\n")
    assert f"## `{CONTIGS / 'full.fa'}`" in partial
    assert f"## `{CONTIGS / 'gap.fa'}`" not in partial


def test_output_truncates_existing_file(datadir: Path, tmp_path: Path) -> None:
    """Given --output pointing at an existing file, When screened, Then the
    run truncates it (v1 overwrite semantics, never append)."""
    out_path = tmp_path / "out.tsv"
    out_path.write_text("stale content that must disappear\n" * 20, encoding="utf-8")
    result = runner.invoke(
        app,
        [
            "screen",
            str(CONTIGS / "full.fa"),
            "--db",
            "tinyamr",
            "--datadir",
            str(datadir),
            "--quiet",
            "--output",
            str(out_path),
        ],
    )
    assert result.exit_code == 0
    assert "stale content" not in out_path.read_text(encoding="utf-8")
    assert out_path.read_text(encoding="utf-8").startswith("#FILE\t")


def test_output_persists_streamed_prefix_on_midbatch_failure(
    datadir: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given --output where file 2 fails mid-batch (monkeypatched blast
    error), When screened, Then file 1's already-streamed rows persist in
    the output file (header + file 1) while the typed envelope exits 5 —
    the documented error interaction for streamed emission."""
    out_path = tmp_path / "out.tsv"
    real = gapit_screening.screen_file

    def fail_on_second(
        query: Path,
        database: Database,
        params: ScreeningParams,
        *,
        dbtype: Literal["nucl", "prot"],
        debug: bool = False,
        merge_fragments: bool = False,
    ) -> Report:
        if query.name == "gap.fa":
            raise InputError(f"boom in {query}", code="INVALID_INPUT", context={})
        return real(
            query, database, params, dbtype=dbtype, debug=debug, merge_fragments=merge_fragments
        )

    monkeypatch.setattr(gapit_screening, "screen_file", fail_on_second)
    result = runner.invoke(
        app,
        [
            "screen",
            str(CONTIGS / "full.fa"),
            str(CONTIGS / "gap.fa"),
            "--db",
            "tinyamr",
            "--datadir",
            str(datadir),
            "--nopath",
            "--output",
            str(out_path),
        ],
    )
    assert result.exit_code == 5
    assert result.stdout == ""
    partial = out_path.read_text(encoding="utf-8")
    assert partial.startswith("#FILE\t")
    assert len(partial.splitlines()) == 2
    assert partial.splitlines()[1].startswith("full.fa\t")


def test_debug_echoes_normalize_and_blast_argv(datadir: Path) -> None:
    """Given --debug, When screening, Then the native normalize step and the
    blastn argv are echoed to stderr (abricate --debug parity; the retired
    any2fasta pipe echo is gone)."""
    result = screen(datadir, "--debug", str(CONTIGS / "full.fa"))
    assert result.exit_code == 0
    normalize_lines = [ln for ln in result.stderr.splitlines() if ln.startswith("gapit: normalize")]
    run_lines = [ln for ln in result.stderr.splitlines() if ln.startswith("gapit: run:")]
    assert normalize_lines == [f"gapit: normalize: {CONTIGS / 'full.fa'} (fasta)"]
    assert len(run_lines) == 1
    assert run_lines[0].startswith("gapit: run: blastn ")
    assert "-perc_identity" in run_lines[0]


def test_debug_stdout_identical_to_plain_run(datadir: Path) -> None:
    """Given the same run with and without --debug, When compared, Then
    stdout is byte-identical (debug is stderr-only)."""
    debug_run = screen(datadir, "--debug", "--nopath", str(CONTIGS / "full.fa"))
    plain_run = screen(datadir, "--nopath", str(CONTIGS / "full.fa"))
    assert debug_run.exit_code == 0
    assert plain_run.exit_code == 0
    assert debug_run.stdout == plain_run.stdout


def test_default_run_emits_no_argv_lines(datadir: Path) -> None:
    """Given a default (no --debug) run, When inspected, Then stderr carries
    no `gapit: run:` argv echo or `gapit: normalize:` lines."""
    result = screen(datadir, str(CONTIGS / "full.fa"))
    assert result.exit_code == 0
    assert "gapit: run:" not in result.stderr
    assert "gapit: normalize:" not in result.stderr


def _screen_with_run_probe_recorder(
    datadir: Path, monkeypatch: pytest.MonkeyPatch, *extra: str
) -> tuple[Result, list[list[str]]]:
    """Run a screen with every subprocess.run argv recorded in order; the real
    binaries still execute (pass-through wrapper over the keyword forms our
    call sites use, including the pipeline's ``input=`` blast invocation)."""
    real_run = subprocess.run
    argvs: list[list[str]] = []

    def counting_run(
        argv: Sequence[str],
        *,
        check: bool = False,
        capture_output: bool = False,
        text: bool = False,
        input: str | None = None,
    ) -> subprocess.CompletedProcess[str]:
        argvs.append(list(argv))
        return real_run(argv, check=check, capture_output=capture_output, text=text, input=input)

    monkeypatch.setattr(subprocess, "run", counting_run)
    return screen(datadir, *extra), argvs


def test_probe_subprocesses_fire_once_across_files(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a 3-file screen, When the per-run probe subprocesses are counted,
    Then `blastn -version` (dependency gate) and `blastdbcmd -info` (dbtype
    resolution) each fire exactly once per run, never once per input file."""
    files = [str(CONTIGS / name) for name in ("full.fa", "gap.fa", "none.fa")]
    result, argvs = _screen_with_run_probe_recorder(datadir, monkeypatch, "--nopath", *files)
    assert result.exit_code == 0
    assert sum(1 for argv in argvs if argv[:2] == ["blastn", "-version"]) == 1
    assert sum(1 for argv in argvs if argv[0] == "blastdbcmd") == 1


def test_probe_subprocesses_fire_once_with_jobs(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the same 3-file screen with --jobs 3, When counted, Then the two
    probes still fire exactly once each — both are hoisted before the pool."""
    files = [str(CONTIGS / name) for name in ("full.fa", "gap.fa", "none.fa")]
    result, argvs = _screen_with_run_probe_recorder(
        datadir, monkeypatch, "--jobs", "3", "--nopath", *files
    )
    assert result.exit_code == 0
    assert sum(1 for argv in argvs if argv[:2] == ["blastn", "-version"]) == 1
    assert sum(1 for argv in argvs if argv[0] == "blastdbcmd") == 1
