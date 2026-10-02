"""--jobs and streaming tests for the positional-FASTQ wildcard path.

The wildcard (``gapit screen -d db *.fastq.gz``) mirrors the contig path's
parallel/streaming contract: ``--jobs N`` screens samples concurrently while
output stays in sample order; tsv (the default) and md stream a static
preamble plus one chunk per completed sample (head-of-line under --jobs);
json stays one document written at the end; ``--output`` persists the
streamed prefix incrementally. The ``--r1``/``--r2`` single-sample path keeps
rejecting ``--jobs`` with the frozen message. Contig-path twins:
test_cli_screen_jobs.py, test_streaming_chunks.py.
"""

import glob
import json
import re
import shutil
import time
from datetime import UTC, datetime
from pathlib import Path

import pytest
from typer.testing import CliRunner, Result

from gapit import screening_reads_positional
from gapit.cli import app
from gapit.db import Database
from gapit.engines import Emit, OutputFormat
from gapit.errors import InputError
from gapit.formats.md import reads_md_chunk, reads_md_preamble
from gapit.formats.reads_tsv import format_reads_tsv, reads_tsv_chunk, reads_tsv_preamble
from gapit.paf import ReadType
from gapit.reads import GeneCoverage, ReadsParams, ReadsReport

READS_DB_DIR = Path(__file__).parent / "data" / "reads_db"
READS = Path(__file__).parent / "data" / "reads"
PINNED_NOW = datetime(2026, 9, 17, 12, 0, 0, tzinfo=UTC)

runner = CliRunner()


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    target = tmp_path / "datadir"
    shutil.copytree(READS_DB_DIR, target)
    return target


@pytest.fixture()
def wildcard_dir(tmp_path: Path) -> Path:
    """Two samples in two marker styles, all .fq.gz (the rightsholder shape)."""
    target = tmp_path / "glob"
    target.mkdir()
    for name in ("s1_1.fq.gz", "s1_2.fq.gz", "s2_R1.fq.gz", "s2_R2.fq.gz"):
        shutil.copy(READS / name, target / name)
    return target


def screen(datadir: Path, *extra: str) -> Result:
    return runner.invoke(app, ["screen", "--db", "tinyreads", "--datadir", str(datadir), *extra])


def _without_timestamps(rendered: str) -> str:
    """created_at differs across runs on a second boundary; comparisons
    normalize it away (the streaming-chunks precedent)."""
    return re.sub(r"created_at: \S+|\"created_at\": \"[^\"]+\"", "<ts>", rendered)


def _use_case(
    datadir: Path,
    wildcard_dir: Path,
    output_format: OutputFormat | None,
    *,
    jobs: int = 1,
    emit: Emit | None = None,
) -> str:
    return screening_reads_positional.run_screen_reads_positional(
        [
            wildcard_dir / name
            for name in ("s1_1.fq.gz", "s1_2.fq.gz", "s2_R1.fq.gz", "s2_R2.fq.gz")
        ],
        "tinyreads",
        datadir,
        None,
        90.0,
        0.0,
        0,
        1,
        output_format,
        True,
        lambda message: None,
        jobs=jobs,
        emit=emit,
    )


# -------------------------------------------------------------- CLI: --jobs --


@pytest.mark.parametrize("flags", [[], ["--format", "json"], ["--format", "md"]])
def test_jobs_2_stdout_identical_to_jobs_1(
    datadir: Path, wildcard_dir: Path, monkeypatch: pytest.MonkeyPatch, flags: list[str]
) -> None:
    """Given the wildcard glob, When screened --jobs 1 vs --jobs 2, Then
    stdout is byte-identical (created_at aside) — samples stay in sample
    order regardless of completion order."""
    monkeypatch.chdir(wildcard_dir)
    files = sorted(glob.glob("*.gz"))
    seq = screen(datadir, "--jobs", "1", *flags, *files)
    par = screen(datadir, "--jobs", "2", *flags, *files)
    assert seq.exit_code == 0
    assert par.exit_code == 0
    assert _without_timestamps(par.stdout) == _without_timestamps(seq.stdout)


def test_jobs_default_matches_explicit_1(
    datadir: Path, wildcard_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the wildcard with no --jobs and with --jobs 1, When compared,
    Then stdout AND stderr are byte-identical (default = sequential path)."""
    monkeypatch.chdir(wildcard_dir)
    files = sorted(glob.glob("*.gz"))
    default = screen(datadir, *files)
    explicit = screen(datadir, "--jobs", "1", *files)
    assert default.exit_code == 0
    assert explicit.exit_code == 0
    assert _without_timestamps(explicit.stdout) == _without_timestamps(default.stdout)
    assert explicit.stderr == default.stderr


def test_jobs_2_stderr_chatter_multiset_matches_jobs_1(
    datadir: Path, wildcard_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the wildcard, When screened --jobs 1 vs --jobs 2, Then stderr
    holds the SAME lines as a multiset — order may interleave under --jobs >
    1 (cpu_count pinned high keeps the oversubscription note out)."""
    monkeypatch.setattr("os.cpu_count", lambda: 64)
    monkeypatch.chdir(wildcard_dir)
    files = sorted(glob.glob("*.gz"))
    seq = screen(datadir, "--jobs", "1", *files)
    par = screen(datadir, "--jobs", "2", *files)
    assert seq.exit_code == 0
    assert par.exit_code == 0
    assert sorted(par.stderr.splitlines()) == sorted(seq.stderr.splitlines())


def test_invalid_jobs_exits_2(
    datadir: Path, wildcard_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given --jobs 0 on the wildcard path, When screened, Then exit 2 with
    a USAGE_ERROR envelope carrying the contig-path message."""
    monkeypatch.chdir(wildcard_dir)
    result = screen(datadir, "--jobs", "0", "s1_1.fq.gz", "s1_2.fq.gz")
    assert result.exit_code == 2
    error = json.loads(result.stderr.splitlines()[-1])
    assert error["code"] == "USAGE_ERROR"
    assert error["message"] == "--jobs must be >= 1: got 0"


def test_jobs_oversubscription_warning(
    datadir: Path, wildcard_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given --jobs 2 on a one-cpu machine (cpu_count pinned), When
    screened, Then the note appears exactly once with the contig-path
    wording, and --quiet suppresses it."""
    monkeypatch.setattr("os.cpu_count", lambda: 1)
    monkeypatch.chdir(wildcard_dir)
    result = screen(datadir, "--jobs", "2", "s1_1.fq.gz", "s1_2.fq.gz")
    assert result.exit_code == 0
    assert "--jobs 2 --threads 1 oversubscribes 1 cpus" in result.stderr
    assert result.stderr.count("oversubscribes") == 1
    quiet = screen(datadir, "--jobs", "2", "--quiet", "s1_1.fq.gz", "s1_2.fq.gz")
    assert quiet.exit_code == 0
    assert quiet.stderr == ""


def test_r1_mode_still_rejects_jobs(
    datadir: Path, wildcard_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given --r1/--r2 (single-sample reads mode) with --jobs 2, When
    screened, Then the frozen rejection still fires (only the wildcard path
    legalized --jobs)."""
    monkeypatch.chdir(wildcard_dir)
    result = screen(datadir, "--r1", "s1_1.fq.gz", "--r2", "s1_2.fq.gz", "--jobs", "2")
    assert result.exit_code == 2
    error = json.loads(result.stderr.splitlines()[-1])
    assert error["message"] == "--jobs is not available in reads mode"


# ------------------------------------------------ --output incremental writes --


def test_output_md_persists_streamed_prefix_stdout_clean(
    datadir: Path, wildcard_dir: Path, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Given the wildcard with --format md -o PATH, When run, Then stdout
    carries no data and the file holds exactly the streamed document
    (preamble first, one section per sample in sample order)."""
    monkeypatch.chdir(wildcard_dir)
    files = sorted(glob.glob("*.gz"))
    target = tmp_path / "out.md"
    to_file = screen(datadir, "--format", "md", "-o", str(target), *files)
    assert to_file.exit_code == 0
    assert to_file.stdout == ""
    to_stdout = screen(datadir, "--format", "md", *files)
    assert _without_timestamps(target.read_text(encoding="utf-8")) == _without_timestamps(
        to_stdout.stdout
    )


def test_output_midbatch_failure_persists_already_written_samples(
    datadir: Path, wildcard_dir: Path, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Given sample 2 failing mid-batch under --format md -o PATH, When run,
    Then the envelope exits nonzero but the file keeps the preamble and
    sample 1's already-streamed section (contig-path semantics); the json
    variant renders nothing, so no artifact appears at all."""
    real = screening_reads_positional.screen_reads

    def fail_s2(
        lanes: list[tuple[Path, Path | None]],
        database: Database,
        *,
        read_type: ReadType,
        min_breadth: float,
        threads: int,
        debug: bool = False,
        min_identity: float = 0.0,
        min_mapq: int = 0,
    ) -> ReadsReport:
        names = {path.name for lane in lanes for path in lane if path is not None}
        if any(name.startswith("s2_") for name in names):
            raise InputError("synthetic sample-2 failure", code="INVALID_READS_FORMAT")
        return real(
            lanes,
            database,
            read_type=read_type,
            min_breadth=min_breadth,
            threads=threads,
            debug=debug,
            min_identity=min_identity,
            min_mapq=min_mapq,
        )

    monkeypatch.setattr(screening_reads_positional, "screen_reads", fail_s2)
    monkeypatch.chdir(wildcard_dir)
    md_target = tmp_path / "out.md"
    md_run = screen(datadir, "--format", "md", "-o", str(md_target), *sorted(glob.glob("*.gz")))
    assert md_run.exit_code == 5
    assert md_run.stdout == ""
    persisted = md_target.read_text(encoding="utf-8")
    assert persisted.startswith("---\nschema: gapit.reads/1\n")
    assert "## `s1`" in persisted
    assert "## `s2`" not in persisted
    json_target = tmp_path / "out.json"
    json_run = screen(
        datadir, "--format", "json", "-o", str(json_target), *sorted(glob.glob("*.gz"))
    )
    assert json_run.exit_code == 5
    assert not json_target.exists()
    tsv_target = tmp_path / "out.tsv"
    tsv_run = screen(datadir, "-o", str(tsv_target), *sorted(glob.glob("*.gz")))
    assert tsv_run.exit_code == 5
    persisted_tsv = tsv_target.read_text(encoding="utf-8")
    assert persisted_tsv.startswith("#SAMPLE\tGENE\tBREADTH%")
    assert "\ns1\t" in persisted_tsv
    assert "\ns2\t" not in persisted_tsv


# ------------------------------------------------- use-case emit contract --


def test_emit_md_chunks_match_buffered_and_arrive_per_sample(
    datadir: Path, wildcard_dir: Path
) -> None:
    """Given a streamed wildcard run (emit=collector) beside a buffered one,
    When compared, Then chunk concatenation equals the returned string and
    the buffered render (created_at aside), the first chunk is the STATIC
    preamble (no run totals), and each sample's section arrives as its own
    chunk in sample order."""
    collected: list[str] = []
    streamed = _use_case(datadir, wildcard_dir, OutputFormat.md, emit=collected.append)
    buffered = _use_case(datadir, wildcard_dir, OutputFormat.md)
    assert streamed == buffered
    assert "".join(collected) == buffered
    assert collected[0].startswith("---\nschema: gapit.reads/1\n")
    assert "files:" not in collected[0] and "genes_found:" not in collected[0]
    assert [chunk.splitlines()[0] for chunk in collected[1:]] == ["## `s1`", "## `s2`"]


def test_emit_tsv_default_chunks_match_buffered_and_arrive_per_sample(
    datadir: Path, wildcard_dir: Path
) -> None:
    """Given a streamed wildcard run with no format (the tsv default), When
    compared beside a buffered explicit-tsv one, Then chunk concatenation
    equals the buffered render, the first chunk is the bare header, and each
    sample's rows arrive as their own chunk in sample order."""
    collected: list[str] = []
    streamed = _use_case(datadir, wildcard_dir, None, emit=collected.append)
    buffered = _use_case(datadir, wildcard_dir, OutputFormat.tsv)
    assert streamed == buffered
    assert "".join(collected) == buffered
    assert collected[0].startswith("#SAMPLE\tGENE\tBREADTH%")
    assert [chunk.splitlines()[0].split("\t")[0] for chunk in collected[1:]] == ["s1"]


def test_emit_json_is_one_chunk_matching_buffered(datadir: Path, wildcard_dir: Path) -> None:
    """Given explicit json (the single-document contract), When streamed
    through emit, Then exactly one chunk equals the buffered render
    (created_at aside) — no partial documents."""
    collected: list[str] = []
    streamed = _use_case(datadir, wildcard_dir, OutputFormat.json, emit=collected.append)
    buffered = _use_case(datadir, wildcard_dir, OutputFormat.json)
    assert _without_timestamps(streamed) == _without_timestamps(buffered)
    assert len(collected) == 1
    assert _without_timestamps(collected[0]) == _without_timestamps(buffered)


def test_sample1_emits_before_slow_sample2_completes(
    datadir: Path, wildcard_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given --jobs 2 with sample 2 artificially slowed, When streaming md,
    Then sample 1's section is emitted while sample 2 is still screening —
    the head-of-line proof (emit fires per completed sample, not at run end)."""
    real = screening_reads_positional.screen_reads
    s2_done: dict[str, float] = {}

    def slow_s2(
        lanes: list[tuple[Path, Path | None]],
        database: Database,
        *,
        read_type: ReadType,
        min_breadth: float,
        threads: int,
        debug: bool = False,
        min_identity: float = 0.0,
        min_mapq: int = 0,
    ) -> ReadsReport:
        report = real(
            lanes,
            database,
            read_type=read_type,
            min_breadth=min_breadth,
            threads=threads,
            debug=debug,
            min_identity=min_identity,
            min_mapq=min_mapq,
        )
        names = {path.name for lane in lanes for path in lane if path is not None}
        if any(name.startswith("s2_") for name in names):
            time.sleep(0.4)
            s2_done["at"] = time.monotonic()
        return report

    monkeypatch.setattr(screening_reads_positional, "screen_reads", slow_s2)
    stamps: list[tuple[str, float]] = []
    _use_case(
        datadir,
        wildcard_dir,
        OutputFormat.md,
        jobs=2,
        emit=lambda chunk: stamps.append((chunk, time.monotonic())),
    )
    s1_at = next(at for chunk, at in stamps if chunk.startswith("## `s1`"))
    assert s1_at < s2_done["at"]


# ------------------------------------------------------ renderer byte-proofs --


def _read_report(sample: str, present: bool) -> ReadsReport:
    """One synthetic sample report — no minimap2 needed for the renderer
    proofs."""
    gene = GeneCoverage(
        database="tinyreads",
        gene="tetX",
        accession="SYN-001",
        function="TETRACYCLINE",
        product="extended resistance determinant tetX",
        tlen=1000,
        breadth_pct=97.7 if present else 40.0,
        mean_depth=2.09,
        reads_mapped=12 if present else 6,
        present=present,
        mean_identity_pct=99.5,
    )
    return ReadsReport(reads=(sample,), genes=(gene,))


def test_reads_md_preamble_is_static_and_chunks_carry_the_schema_variants() -> None:
    """Given synthetic reports, When split into preamble + per-sample chunks,
    Then the preamble holds only pre-sample-1 facts (reads/2 adds the two
    filter lines) and each chunk is exactly its sample's section (reads/2
    adds the Identity% column) — the md_report_preamble/chunk split,
    reads edition."""
    params = ReadsParams(db="tinyreads")
    report = _read_report("s1", present=True)
    preamble = reads_md_preamble(params, now=PINNED_NOW)
    assert preamble.startswith("---\nschema: gapit.reads/1\n")
    assert preamble.endswith("# gapit read screening report\n\n")
    assert "files:" not in preamble and "genes_found:" not in preamble
    reads2_preamble = reads_md_preamble(params, now=PINNED_NOW, reads2=True)
    assert "schema: gapit.reads/2" in reads2_preamble
    assert "min_identity: 0.0" in reads2_preamble and "min_mapq: 0" in reads2_preamble
    chunk = reads_md_chunk(report)
    assert chunk.splitlines()[0] == "## `s1`"
    assert "| tetX |" in chunk
    assert "99.50" not in chunk  # reads/1 carries no Identity% column
    reads2_chunk = reads_md_chunk(report, reads2=True)
    assert reads2_chunk.splitlines()[0] == "## `s1`"
    assert "99.50" in reads2_chunk


def test_reads_tsv_preamble_and_chunks_match_format_reads_tsv() -> None:
    """Given synthetic reports, When split into preamble + per-sample chunks,
    Then the concatenation is exactly format_reads_tsv, reads/2 adds the
    Identity% column (header and cells), and csv swaps the separator — the
    reads_md_preamble/chunk split, tsv edition."""
    reports = [_read_report("s1", present=True), _read_report("s2", present=False)]
    composed = reads_tsv_preamble() + "".join(reads_tsv_chunk(r) for r in reports)
    assert composed == format_reads_tsv(reports)
    assert composed.splitlines()[0] == (
        "#SAMPLE\tGENE\tBREADTH%\tDEPTH\tREADS\tPRESENT\tDATABASE\tACCESSION\tPRODUCT"
    )
    assert "s1\ttetX\t97.70\t2.09\t12\tyes\t" in composed
    assert "s2" not in composed  # absent calls hidden by default
    all_composed = reads_tsv_preamble() + "".join(
        reads_tsv_chunk(r, all_genes=True) for r in reports
    )
    assert "s2\ttetX\t40.00\t2.09\t6\tno\t" in all_composed
    assert "\t99.50\t" not in composed
    reads2_composed = reads_tsv_preamble(reads2=True) + "".join(
        reads_tsv_chunk(r, reads2=True) for r in reports
    )
    assert reads2_composed == format_reads_tsv(reports, reads2=True)
    assert reads2_composed.splitlines()[0].endswith("\tIDENTITY%")
    assert "\t99.50" in reads2_composed
    assert format_reads_tsv(reports, csv=True).splitlines()[0].startswith("#SAMPLE,GENE,")
