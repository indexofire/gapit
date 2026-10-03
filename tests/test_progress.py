"""Tests for the shared stderr progress bar (gapit.progress).

The bar is a stderr-side UX layer only: stdout stays byte-identical with
the bar on vs off, --quiet and non-TTY stderr keep today's legacy per-file
notes (no-op shim), and an interactive stderr gets the gmlst column set
with the per-file notes suppressed.
"""

import io
import re
import shutil
import sys
from dataclasses import replace
from pathlib import Path

import pytest
from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)

import gapit.progress
from gapit.db import make_blast_db
from gapit.db_ops import ProviderReceipt, perform_fetch
from gapit.progress import IdleBatch, make_progress, screen_progress
from gapit.providers import REGISTRY
from gapit.screening import OutputFormat, run_screen
from gapit.screening_reads_positional import run_screen_reads_positional

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"
CONTIGS = Path(__file__).parent / "data" / "contigs"
READS_DB_DIR = Path(__file__).parent / "data" / "reads_db"
READS = Path(__file__).parent / "data" / "reads"
RESFINDER_ZIP = Path(__file__).parent / "data" / "bundled_differential" / "resfinder" / "HEAD.zip"


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    target = tmp_path / "datadir"
    shutil.copytree(FIXTURE_DB_DIR, target)
    make_blast_db(target / "tinyamr" / "sequences", "tinyamr")
    return target


@pytest.fixture()
def reads_datadir(tmp_path: Path) -> Path:
    target = tmp_path / "datadir"
    shutil.copytree(READS_DB_DIR, target)
    return target


def tty_console() -> Console:
    """A recording console that fakes an interactive stderr."""
    return Console(file=io.StringIO(), record=True, force_terminal=True, width=100)


def pipe_console() -> Console:
    """A recording console that fakes piped (non-TTY) stderr."""
    return Console(file=io.StringIO(), record=True, force_terminal=False, width=100)


def test_make_progress_uses_gmlst_columns_on_stderr() -> None:
    """Given make_progress, When constructed, Then the gmlst column set
    (spinner, description, bar, M/N, percent, elapsed) rides the
    stderr-bound console (stdout purity)."""
    progress = make_progress()
    assert [type(column) for column in progress.columns] == [
        SpinnerColumn,
        TextColumn,
        BarColumn,
        MofNCompleteColumn,
        TaskProgressColumn,
        TimeElapsedColumn,
    ]
    assert progress.console.file is sys.stderr


def test_quiet_yields_noop_shim(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --quiet, When the batch context opens, Then the shim is yielded
    with no task id (call sites keep printing the legacy notes)."""
    monkeypatch.setattr(gapit.progress, "STATUS_CONSOLE", tty_console())
    with screen_progress(["a.fa", "b.fa"], quiet=True) as (bar, task_id):
        assert task_id is None
        assert isinstance(bar, IdleBatch)
        bar.describe("a.fa")
        bar.advance()


def test_non_tty_stderr_yields_noop_shim(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given stderr is not a TTY (piped/CI), When quiet is off, Then the bar
    stays disabled — logs keep today's per-file notes, never bar frames."""
    sink = pipe_console()
    monkeypatch.setattr(gapit.progress, "STATUS_CONSOLE", sink)
    with screen_progress(["a.fa", "b.fa"], quiet=False) as (bar, task_id):
        assert task_id is None
        assert isinstance(bar, IdleBatch)
        bar.describe("a.fa")
        bar.advance()
    assert sink.export_text(clear=False).strip() == ""


def test_tty_stderr_activates_the_bar(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given stderr is a TTY and quiet is off, When the batch advances, Then
    a real task is yielded and the captured frames show the gmlst columns
    (description, M/N, percent, elapsed)."""
    sink = tty_console()
    monkeypatch.setattr(gapit.progress, "STATUS_CONSOLE", sink)
    with screen_progress(["a.fa", "b.fa"], quiet=False) as (bar, task_id):
        assert task_id is not None
        bar.describe("a.fa")
        bar.advance()
        bar.describe("b.fa")
        bar.advance()
    text = sink.export_text(clear=False)
    assert "Screening" in text
    assert "2/2" in text and "100%" in text and "0:00:00" in text


def test_bar_renders_current_description_and_counts(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given the bar mid-batch, When a file starts and one completes, Then a
    forced refresh renders the current file's basename and 1/2."""
    sink = tty_console()
    monkeypatch.setattr(gapit.progress, "STATUS_CONSOLE", sink)
    progress = make_progress()
    with progress:
        task_id = progress.add_task("Screening", total=2)
        progress.update(task_id, description="full.fa")
        progress.advance(task_id)
        progress.refresh()
    text = sink.export_text(clear=False)
    assert "full.fa" in text and "1/2" in text


def screen_once(datadir: Path, files: list[Path], output_format: OutputFormat) -> str:
    """One blastn contig batch through the use-case (the CLI's callee)."""
    return run_screen(
        files,
        "tinyamr",
        datadir,
        80.0,
        80.0,
        1,
        1,
        fofn=None,
        quiet=False,
        noheader=False,
        nopath=True,
        debug=False,
        output_format=output_format,
    )


def test_screen_bar_keeps_stdout_byte_identical_and_suppresses_notes(
    datadir: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Given two contigs, When screened with the bar active vs inactive, Then
    stdout is byte-identical, notes print ONLY when the bar is off, and the
    captured bar carries the batch frames."""
    files = [CONTIGS / "full.fa", CONTIGS / "gap.fa"]
    plain = screen_once(datadir, files, OutputFormat.tsv)
    plain_err = capsys.readouterr().err

    sink = tty_console()
    monkeypatch.setattr(gapit.progress, "STATUS_CONSOLE", sink)
    barred = screen_once(datadir, files, OutputFormat.tsv)
    barred_err = capsys.readouterr().err

    assert barred == plain
    assert "Processing:" in plain_err and "Found" in plain_err
    assert "Processing:" not in barred_err and "Found" not in barred_err
    bar_text = sink.export_text(clear=False)
    assert "Screening" in bar_text and "2/2" in bar_text


def test_screen_md_stream_byte_identical_with_bar_on(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given --format md (streaming frontmatter + per-file sections), When
    screened with the bar active vs inactive, Then the documents are
    identical modulo created_at (the frontmatter timestamp)."""
    files = [CONTIGS / "full.fa", CONTIGS / "gap.fa"]
    plain = screen_once(datadir, files, OutputFormat.md)
    monkeypatch.setattr(gapit.progress, "STATUS_CONSOLE", tty_console())
    barred = screen_once(datadir, files, OutputFormat.md)
    stamp = re.compile(r"created_at: \d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z")
    assert stamp.sub("created_at: TS", barred) == stamp.sub("created_at: TS", plain)


def reads_once(datadir: Path, files: list[Path]) -> str:
    """One positional-FASTQ batch through the use-case (the CLI's callee)."""

    def warn(message: str) -> None:
        raise AssertionError(f"unexpected pairing warning: {message}")

    return run_screen_reads_positional(
        files,
        "tinyreads",
        datadir,
        read_type=None,
        min_breadth=90.0,
        min_identity=0.0,
        min_mapq=0,
        threads=1,
        output_format=OutputFormat.tsv,
        quiet=False,
        warn=warn,
    )


def test_reads_bar_keeps_streaming_header_first_and_suppresses_notes(
    reads_datadir: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Given the four-file two-sample wildcard shape, When screened with the
    bar active vs inactive, Then stdout is byte-identical and still leads
    with the lazy #SAMPLE header, and the per-sample notes print ONLY when
    the bar is off."""
    files = [
        READS / "s1_1.fq.gz",
        READS / "s1_2.fq.gz",
        READS / "s2_R1.fq.gz",
        READS / "s2_R2.fq.gz",
    ]
    plain = reads_once(reads_datadir, files)
    plain_err = capsys.readouterr().err

    sink = tty_console()
    monkeypatch.setattr(gapit.progress, "STATUS_CONSOLE", sink)
    barred = reads_once(reads_datadir, files)
    barred_err = capsys.readouterr().err

    assert barred == plain
    assert barred.startswith("#SAMPLE\t")
    assert "Screening sample s1 reads:" in plain_err
    assert "Screening sample" not in barred_err and "Detected" not in barred_err
    bar_text = sink.export_text(clear=False)
    assert "Screening" in bar_text and "2/2" in bar_text


def test_fetch_bar_advances_per_db_and_keeps_receipts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given `db fetch all` over two offline (file://) providers, When
    fetched with the bar active vs quiet, Then the receipts carry the same
    (db, records, dbtype) and the captured bar advances 2/2."""
    provider = REGISTRY["resfinder"]
    local = replace(provider, source_urls=(RESFINDER_ZIP.as_uri(),))
    monkeypatch.setattr("gapit.db_ops.REGISTRY", {"resfinder": local, "resfinder2": local})
    monkeypatch.setattr("gapit.db_ops.DEFAULT_DBS", ("resfinder", "resfinder2"))

    quiet_dd = tmp_path / "quiet-dd"
    quiet_dd.mkdir()
    quiet_receipts = list(perform_fetch("all", quiet_dd, quiet=True))

    bar_dd = tmp_path / "bar-dd"
    bar_dd.mkdir()
    sink = tty_console()
    monkeypatch.setattr(gapit.progress, "STATUS_CONSOLE", sink)
    bar_receipts = list(perform_fetch("all", bar_dd, quiet=False))

    def shape(receipts: list[ProviderReceipt]) -> list[tuple[str, int, str]]:
        return [(receipt.db, receipt.records, str(receipt.dbtype)) for receipt in receipts]

    assert shape(bar_receipts) == shape(quiet_receipts)
    assert [receipt.db for receipt in bar_receipts] == ["resfinder", "resfinder2"]
    bar_text = sink.export_text(clear=False)
    assert "Fetching" in bar_text and "2/2" in bar_text
