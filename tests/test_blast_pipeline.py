"""Pipeline stderr-drain tests for blast._pipeline (#6): any2fasta stderr is
drained by a background thread — an undrained stderr pipe fills at ~64 KB and
would deadlock the `any2fasta | blastn` pipeline (the classic Popen deadlock;
watchdog pattern from tests/test_reads_streaming.py). Fakes are executable
python scripts named ``any2fasta``/``blastn`` on a monkeypatched PATH.
"""

import os
import threading
from pathlib import Path

import pytest

from gapit.blast import screen_file
from gapit.db import Database
from gapit.errors import InputError
from gapit.report import ScreeningParams

# > 64 KB so the child's stderr write blocks unless something drains it.
STDERR_SPEW = "spew-line-that-never-ends\n" * 3000

FASTA_OUT = ">contig1\nACGTACGTAC\n"


def write_fake(
    bin_dir: Path, name: str, *, stderr_text: str, stdout_text: str, exit_code: int
) -> None:
    """An executable fake that spews stderr, then writes stdout, then exits
    with exit_code. The stderr write precedes stdout so a run that does not
    drain stderr concurrently can never produce output."""
    script = (
        "#!/usr/bin/env python3\n"
        "import sys\n"
        f"sys.stderr.write({stderr_text!r})\n"
        "sys.stderr.flush()\n"
        f"sys.stdout.write({stdout_text!r})\n"
        "sys.stdout.flush()\n"
        f"sys.exit({exit_code})\n"
    )
    binary = bin_dir / name
    binary.write_text(script, encoding="utf-8")
    binary.chmod(0o755)


def write_fake_blastn(bin_dir: Path) -> None:
    """An executable fake blastn: consume stdin fully (the pipeline contract),
    emit no rows, exit 0."""
    script = "#!/usr/bin/env python3\nimport sys\nsys.stdin.read()\nsys.exit(0)\n"
    binary = bin_dir / "blastn"
    binary.write_text(script, encoding="utf-8")
    binary.chmod(0o755)


@pytest.fixture()
def fakedb(tmp_path: Path) -> Database:
    """A Database whose paths are never read (the fake binaries ignore argv)."""
    return Database(name="fakedb", path=tmp_path, sequences_path=tmp_path / "sequences")


def fake_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")
    return bin_dir


def call_with_watchdog(query: Path, database: Database) -> tuple[int, InputError | None]:
    """Screen one file in a daemon thread; a pipeline deadlock becomes a clean
    assertion failure instead of a hung suite, and a typed error raised inside
    the thread is re-raised on the main thread."""
    params = ScreeningParams(db="fakedb", minid=80.0, mincov=80.0, threads=1)
    hits: list[int] = []
    failure: list[Exception] = []

    def call() -> None:
        try:
            report = screen_file(query, database, params, dbtype="nucl")
            hits.append(len(report.hits))
        except Exception as exc:  # transport to the joining thread
            failure.append(exc)

    worker = threading.Thread(target=call, daemon=True)
    worker.start()
    worker.join(timeout=30)
    assert not worker.is_alive(), "pipeline deadlocked on an undrained any2fasta stderr pipe"
    if failure:
        assert isinstance(failure[0], InputError)
        return hits[0] if hits else -1, failure[0]
    return hits[0], None


def test_any2fasta_stderr_spew_does_not_deadlock(
    tmp_path: Path, fakedb: Database, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given an any2fasta that writes >64 KB of stderr BEFORE its FASTA stdout
    and a blastn that consumes all of stdin, When screened, Then the pipeline
    completes with zero hits (fake blastn emits no rows)."""
    bin_dir = fake_path(monkeypatch, tmp_path)
    write_fake(bin_dir, "any2fasta", stderr_text=STDERR_SPEW, stdout_text=FASTA_OUT, exit_code=0)
    write_fake_blastn(bin_dir)
    query = tmp_path / "contigs.fa"
    query.write_text(FASTA_OUT, encoding="utf-8")
    n_hits, error = call_with_watchdog(query, fakedb)
    assert error is None
    assert n_hits == 0


def test_any2fasta_failure_after_spew_carries_stderr_text(
    tmp_path: Path, fakedb: Database, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given an any2fasta that spews stderr and exits 1 while blastn exits 0,
    When screened, Then INVALID_INPUT carries the drained stderr text (the
    background drain feeds the error message, not just the deadlock fix)."""
    bin_dir = fake_path(monkeypatch, tmp_path)
    write_fake(bin_dir, "any2fasta", stderr_text=STDERR_SPEW, stdout_text=FASTA_OUT, exit_code=1)
    write_fake_blastn(bin_dir)
    query = tmp_path / "contigs.fa"
    query.write_text(FASTA_OUT, encoding="utf-8")
    _, error = call_with_watchdog(query, fakedb)
    assert error is not None
    assert error.code == "INVALID_INPUT"
    assert "spew-line-that-never-ends" in str(error)
    assert error.context["file"] == str(query)
