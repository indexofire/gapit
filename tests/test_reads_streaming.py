"""Streaming tests for run_minimap2 (#2): stdout is consumed line by line
while a background thread drains stderr — an undrained stderr pipe fills at
~64 KB and would block the child mid-run (the classic Popen deadlock).
Fakes are executable python scripts named ``minimap2`` on a monkeypatched
PATH (argv stays a list; shell is never involved).
"""

import os
import subprocess
import threading
from io import StringIO
from pathlib import Path

import pytest

from gapit.db import Database
from gapit.errors import DependencyError, GapitError
from gapit.paf import PafRecord
from gapit.reads import run_minimap2

# > 64 KB so the child's stderr write blocks unless something drains it.
STDERR_SPEW = "spew-line-that-never-ends\n" * 3000

PAF_ROWS = (
    "r1\t100\t0\t50\t+\tdb~~~geneA~~~ACC~~~RES\t100\t10\t60\t50\t50\t60\ttp:A:P\n"
    "r2\t100\t0\t50\t+\tdb~~~geneA~~~ACC~~~RES\t100\t70\t90\t20\t20\t60\n"
)


def write_fake_minimap2(
    bin_dir: Path, *, stderr_text: str, stdout_text: str, exit_code: int
) -> None:
    """An executable ``minimap2`` that spews stderr, then writes stdout, then
    exits with exit_code. The stderr write precedes stdout so a run that does
    not drain stderr concurrently can never reach the PAF rows."""
    script = (
        "#!/usr/bin/env python3\n"
        "import sys\n"
        f"sys.stderr.write({stderr_text!r})\n"
        "sys.stderr.flush()\n"
        f"sys.stdout.write({stdout_text!r})\n"
        "sys.stdout.flush()\n"
        f"sys.exit({exit_code})\n"
    )
    binary = bin_dir / "minimap2"
    binary.write_text(script, encoding="utf-8")
    binary.chmod(0o755)


@pytest.fixture()
def fakedb(tmp_path: Path) -> Database:
    """A Database whose paths are never read (the fake binary ignores argv)."""
    return Database(name="fakedb", path=tmp_path, sequences_path=tmp_path / "sequences")


def fake_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")
    return bin_dir


def call_with_watchdog(database: Database, tmp_path: Path) -> list[PafRecord]:
    """Run one lane in a daemon thread; a streaming deadlock becomes a
    clean assertion failure instead of a hung suite, and a typed error
    raised inside the thread is re-raised on the main thread."""
    reads = tmp_path / "reads.fq"
    reads.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")
    outcome: list[PafRecord] = []
    failure: list[Exception] = []

    def call() -> None:
        try:
            outcome.extend(run_minimap2([(reads, None)], database, read_type="sr", threads=1))
        except Exception as exc:  # transport to the joining thread
            failure.append(exc)

    worker = threading.Thread(target=call, daemon=True)
    worker.start()
    worker.join(timeout=30)
    assert not worker.is_alive(), "run_minimap2 deadlocked on an undrained stderr pipe"
    if failure:
        raise failure[0]
    return outcome


def test_stderr_spew_mid_run_does_not_deadlock(
    tmp_path: Path, fakedb: Database, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a minimap2 that writes >64 KB of stderr BEFORE its PAF stdout,
    When run, Then the invocation completes, both rows parse, and the tp:A:P
    / no-tp primary logic is preserved."""
    bin_dir = fake_path(monkeypatch, tmp_path)
    write_fake_minimap2(bin_dir, stderr_text=STDERR_SPEW, stdout_text=PAF_ROWS, exit_code=0)
    rows = call_with_watchdog(fakedb, tmp_path)
    assert [row.qname for row in rows] == ["r1", "r2"]
    assert rows[0].is_primary is True
    assert rows[1].is_primary is True
    assert rows[0].tname == "db~~~geneA~~~ACC~~~RES"


def test_nonzero_exit_after_spew_raises_typed_error(
    tmp_path: Path, fakedb: Database, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a minimap2 that spews stderr and exits 1, When run, Then
    MINIMAP2_FAILED carries the captured stderr text."""
    bin_dir = fake_path(monkeypatch, tmp_path)
    write_fake_minimap2(bin_dir, stderr_text=STDERR_SPEW, stdout_text="", exit_code=1)
    with pytest.raises(GapitError) as excinfo:
        call_with_watchdog(fakedb, tmp_path)
    assert excinfo.value.code == "MINIMAP2_FAILED"
    assert "spew-line-that-never-ends" in str(excinfo.value)
    assert excinfo.value.context["binary"] == "minimap2"


def test_missing_binary_raises_dependency_error(
    tmp_path: Path, fakedb: Database, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a PATH without minimap2, When run, Then typed DependencyError
    (exit 3 class) with the binary in context."""
    empty_bin = tmp_path / "emptybin"
    empty_bin.mkdir()
    monkeypatch.setenv("PATH", str(empty_bin))
    reads = tmp_path / "reads.fq"
    reads.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")
    with pytest.raises(DependencyError) as excinfo:
        run_minimap2([(reads, None)], fakedb, read_type="sr", threads=1)
    assert excinfo.value.code == "MISSING_DEPENDENCY"
    assert excinfo.value.context["binary"] == "minimap2"


class _FakePopen:
    """The Popen surface _stream_minimap2 touches: stdout line iteration,
    stderr read, wait, returncode."""

    def __init__(self) -> None:
        self.stdout = StringIO(PAF_ROWS)
        self.stderr = StringIO()
        self.returncode = 0

    def wait(self) -> int:
        return 0


def test_input_paths_are_absolute_and_stdin_is_devnull(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a query file literally named '-d' paired with a mate (the
    index-dump injection shape) and a relative db sequences path, When run,
    Then every path argument in argv is absolute — '/…/-d' cannot parse as
    the -d option — and the child's stdin is /dev/null, never an inherited
    stdin (the MCP server's JSON-RPC stream)."""
    database = Database(name="fakedb", path=Path("fakedb"), sequences_path=Path("sequences"))
    reads = tmp_path / "-d"
    reads.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")
    mate = tmp_path / "mate.fq"
    mate.write_text("@r2\nACGT\n+\nIIII\n", encoding="utf-8")
    calls: list[tuple[list[str], dict[str, object]]] = []

    def fake_popen(argv: list[str], **kwargs: object) -> _FakePopen:
        calls.append((argv, kwargs))
        return _FakePopen()

    monkeypatch.setattr(subprocess, "Popen", fake_popen)
    rows = run_minimap2([(reads, mate)], database, read_type="sr", threads=1)
    assert [row.qname for row in rows] == ["r1", "r2"]
    ((argv, kwargs),) = calls
    assert argv[:5] == ["minimap2", "-x", "sr", "-t", "1"]
    assert kwargs["stdin"] == subprocess.DEVNULL
    db_arg, r1_arg, r2_arg = argv[5:8]
    assert Path(db_arg) == Path("sequences").absolute()
    assert Path(r1_arg) == reads.absolute()
    assert Path(r2_arg) == mate.absolute()
    assert not any(argument.startswith("-") for argument in argv[5:])
