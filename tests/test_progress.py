"""Tests for the shared stderr progress bar (gapit.progress).

The bar is a stderr-side UX layer only: stdout stays byte-identical with
the bar on vs off, --quiet and non-TTY stderr keep today's legacy per-file
notes (no-op shim), and an interactive stderr gets the gmlst column set
with the per-file notes suppressed.
"""

import io
import os
import re
import shutil
import sys
from collections.abc import Iterable
from pathlib import Path
from types import SimpleNamespace

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
from gapit.fasta import iter_fasta
from gapit.progress import IdleBatch, make_progress, screen_progress
from gapit.providers.common import Provider
from gapit.records import Record
from gapit.screening import OutputFormat, run_screen
from gapit.screening_reads_positional import run_screen_reads_positional

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"
CONTIGS = Path(__file__).parent / "data" / "contigs"
READS_DB_DIR = Path(__file__).parent / "data" / "reads_db"
READS = Path(__file__).parent / "data" / "reads"


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


def _fixed_terminal_size(_fallback: tuple[int, int]) -> os.terminal_size:
    return os.terminal_size((60, 20))


def tty_console() -> Console:
    """A recording console that fakes an interactive stderr."""
    return Console(file=io.StringIO(), record=True, force_terminal=True, width=100)


def pipe_console() -> Console:
    """A recording console that fakes piped (non-TTY) stderr."""
    return Console(file=io.StringIO(), record=True, force_terminal=False, width=100)


def synthetic_provider(tmp_path: Path, name: str) -> Provider:
    """A tiny offline provider whose source is a local FASTA fixture."""
    source = tmp_path / f"{name}.fa"
    source.write_text(">tetA\n" + "ACGTAG" * 10 + "\n", encoding="utf-8")

    def transform(workdir: Path) -> Iterable[Record]:
        for record in iter_fasta(workdir / source.name):
            yield Record(db=name, gene=record.id, sequence=record.sequence)

    return Provider(
        name=name,
        description=f"synthetic {name} provider for progress tests",
        vendor="Synthetica",
        source_urls=(source.as_uri(),),
        dbtype="nucl",
        transform=transform,
    )


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
    """Given `db fetch -d all` over two offline (file://) providers, When
    fetched with the bar active vs quiet, Then the receipts carry the same
    (db, records, dbtype) and the captured bar advances 2/2."""
    first = synthetic_provider(tmp_path, "synone")
    second = synthetic_provider(tmp_path, "syntwo")
    monkeypatch.setattr("gapit.db_ops.REGISTRY", {first.name: first, second.name: second})
    monkeypatch.setattr("gapit.db_ops.DEFAULT_DBS", (first.name, second.name))

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
    assert [receipt.db for receipt in bar_receipts] == ["synone", "syntwo"]
    bar_text = sink.export_text(clear=False)
    assert "Fetching" in bar_text and "2/2" in bar_text


def test_live_batch_suspend_stops_and_restarts_the_bar(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Given a live bar, When stdout data is emitted under suspend(), Then
    the bar is stopped, the write happens on a clean line, and the bar is
    restarted and refreshed."""
    events: list[str] = []

    class FakeProgress:
        console = SimpleNamespace(is_terminal=True)

        def add_task(self, description: str, *, total: int) -> int:
            events.append(f"add:{description}")
            return 1

        def update(self, task_id: int, *, description: str) -> None:
            events.append(f"describe:{description}")

        def advance(self, task_id: int) -> None:
            events.append("advance")

        def stop(self) -> None:
            events.append("stop")

        def start(self) -> None:
            events.append("start")

        def refresh(self) -> None:
            events.append("refresh")

        def __enter__(self) -> "FakeProgress":
            return self

        def __exit__(self, *args: object) -> None:
            return None

    monkeypatch.setattr(gapit.progress, "make_progress", FakeProgress)
    with screen_progress(["a", "b"], quiet=False) as (bar, task_id):
        assert task_id == 1
        with bar.suspend():
            events.append("data")
    assert events == ["add:Screening", "stop", "data", "start"]


def test_fetch_emit_callback_runs_inside_the_bar(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given perform_fetch(emit=...), When receipt lines are emitted, Then the
    same receipts are yielded and the callback sees every database."""
    first = synthetic_provider(tmp_path, "synone")
    second = synthetic_provider(tmp_path, "syntwo")
    monkeypatch.setattr("gapit.db_ops.REGISTRY", {first.name: first, second.name: second})
    monkeypatch.setattr("gapit.db_ops.DEFAULT_DBS", (first.name, second.name))
    datadir = tmp_path / "dd"
    datadir.mkdir()
    emitted: list[ProviderReceipt] = []
    received = list(perform_fetch("all", datadir, quiet=True, emit=emitted.append))
    assert [r.db for r in emitted] == ["synone", "syntwo"]
    assert emitted == received


def test_live_batch_write_stops_the_bar_for_stdout(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Given a live bar, When data is written through write(), Then the bar
    stops, stdout receives the bytes on a clean line, and the bar restarts."""
    events: list[str] = []

    class FakeProgress:
        console = SimpleNamespace(is_terminal=True)

        def add_task(self, description: str, *, total: int) -> int:
            return 1

        def update(self, task_id: int, *, description: str) -> None:
            return None

        def advance(self, task_id: int) -> None:
            return None

        def stop(self) -> None:
            events.append("stop")

        def start(self) -> None:
            events.append("start")

        def __enter__(self) -> "FakeProgress":
            return self

        def __exit__(self, *args: object) -> None:
            return None

    monkeypatch.setattr(gapit.progress, "make_progress", FakeProgress)
    with screen_progress(["a", "b"], quiet=False) as (bar, _task_id):
        bar.write("row\n")
    assert events == ["stop", "start"]
    assert capsys.readouterr().out == "row\n"


def test_terminal_rule_is_tty_only(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given a pipe, When terminal_rule is requested, Then it is empty; on a
    TTY it is one full-width line, so redirected output stays unchanged."""
    monkeypatch.setattr(
        gapit.progress, "sys", SimpleNamespace(stdout=SimpleNamespace(isatty=lambda: False))
    )
    assert gapit.progress.terminal_rule() == ""

    monkeypatch.setattr(
        gapit.progress, "sys", SimpleNamespace(stdout=SimpleNamespace(isatty=lambda: True))
    )
    monkeypatch.setattr(
        gapit.progress.shutil,
        "get_terminal_size",
        _fixed_terminal_size,
    )
    assert gapit.progress.terminal_rule() == "─" * 60 + "\n"


def test_screen_separates_samples_only_on_interactive_tsv(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given two contigs on an interactive TSV stream, When screened, Then a
    horizontal rule follows the header, separates the file blocks, and closes
    the table — but never appears in CSV or piped output."""
    monkeypatch.setattr(
        gapit.progress, "sys", SimpleNamespace(stdout=SimpleNamespace(isatty=lambda: True))
    )
    monkeypatch.setattr(
        gapit.progress.shutil,
        "get_terminal_size",
        _fixed_terminal_size,
    )
    rule = "─" * 60 + "\n"
    files = [CONTIGS / "full.fa", CONTIGS / "gap.fa"]

    joined = run_screen(
        files,
        "tinyamr",
        datadir,
        80.0,
        80.0,
        1,
        1,
        fofn=None,
        quiet=True,
        noheader=False,
        nopath=True,
        debug=False,
        output_format=OutputFormat.tsv,
        stdout_emit=True,
    )
    assert joined.count(rule) == 3
    header_rule = joined.index(rule)
    between_rule = joined.index(rule, header_rule + 1)
    trailing_rule = joined.rindex(rule)
    assert joined.startswith("#FILE")
    assert header_rule < joined.index("full.fa") < between_rule < joined.index("gap.fa")
    assert trailing_rule > joined.index("gap.fa")
    assert joined.endswith(rule)

    csv_out = run_screen(
        files,
        "tinyamr",
        datadir,
        80.0,
        80.0,
        1,
        1,
        fofn=None,
        quiet=True,
        noheader=False,
        nopath=True,
        debug=False,
        output_format=OutputFormat.csv,
        stdout_emit=True,
    )
    assert rule not in csv_out

    piped = run_screen(
        files,
        "tinyamr",
        datadir,
        80.0,
        80.0,
        1,
        1,
        fofn=None,
        quiet=True,
        noheader=False,
        nopath=True,
        debug=False,
        output_format=OutputFormat.tsv,
    )
    assert rule not in piped


def test_screen_terminal_rules_for_a_single_sample_and_noheader(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given one interactive TSV sample, Then only the header rule and the
    trailing rule appear; --noheader drops the header rule but keeps the
    trailing one."""
    monkeypatch.setattr(
        gapit.progress, "sys", SimpleNamespace(stdout=SimpleNamespace(isatty=lambda: True))
    )
    monkeypatch.setattr(
        gapit.progress.shutil,
        "get_terminal_size",
        _fixed_terminal_size,
    )
    rule = "─" * 60 + "\n"

    joined = run_screen(
        [CONTIGS / "full.fa"],
        "tinyamr",
        datadir,
        80.0,
        80.0,
        1,
        1,
        fofn=None,
        quiet=True,
        noheader=False,
        nopath=True,
        debug=False,
        output_format=OutputFormat.tsv,
        stdout_emit=True,
    )
    assert joined.count(rule) == 2
    assert joined.startswith("#FILE")
    assert joined.endswith(rule)

    noheader = run_screen(
        [CONTIGS / "full.fa"],
        "tinyamr",
        datadir,
        80.0,
        80.0,
        1,
        1,
        fofn=None,
        quiet=True,
        noheader=True,
        nopath=True,
        debug=False,
        output_format=OutputFormat.tsv,
        stdout_emit=True,
    )
    assert noheader.count(rule) == 1
    assert not noheader.startswith(rule)
    assert noheader.endswith(rule)
