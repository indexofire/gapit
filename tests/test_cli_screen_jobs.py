"""Integration tests: --jobs parallel screening (gapit extension) — input-order
stdout, threads determinism, and validation, over the real pipeline on tinyamr."""

import json
import shutil
from pathlib import Path

import pytest
from typer.testing import CliRunner, Result

from gapit.cli import app
from gapit.db import make_blast_db

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"
CONTIGS = Path(__file__).parent / "data" / "contigs"
# Six inputs with duplicate paths (dedup only exists in summary mode).
SIX_FILES = [
    CONTIGS / name for name in ("full.fa", "sort.fa", "gap.fa", "full.fa", "none.fa", "sort.fa")
]

runner = CliRunner()


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


def test_jobs_4_stdout_identical_to_jobs_1(datadir: Path) -> None:
    """Given six inputs (duplicates included), When screened --jobs 1 vs
    --jobs 4, Then stdout is byte-identical (SPEC §4: input order always)."""
    files = [str(path) for path in SIX_FILES]
    seq = screen(datadir, "--jobs", "1", "--nopath", *files)
    par = screen(datadir, "--jobs", "4", "--nopath", *files)
    assert seq.exit_code == 0
    assert par.exit_code == 0
    assert par.stdout == seq.stdout


def test_jobs_default_matches_explicit_1(datadir: Path) -> None:
    """Given the same inputs with no --jobs and with --jobs 1, When compared,
    Then stdout AND stderr are byte-identical (default = sequential path)."""
    files = [str(path) for path in SIX_FILES]
    default = screen(datadir, "--nopath", *files)
    explicit = screen(datadir, "--jobs", "1", "--nopath", *files)
    assert default.exit_code == 0
    assert explicit.exit_code == 0
    assert explicit.stdout == default.stdout
    assert explicit.stderr == default.stderr


def test_jobs_4_stderr_chatter_multiset_matches_jobs_1(datadir: Path) -> None:
    """Given six inputs, When screened --jobs 1 vs --jobs 4, Then stderr holds
    the SAME lines as a multiset — order may interleave under --jobs > 1, so
    lines are compared sorted, never in sequence."""
    files = [str(path) for path in SIX_FILES]
    seq = screen(datadir, "--jobs", "1", *files)
    par = screen(datadir, "--jobs", "4", *files)
    assert seq.exit_code == 0
    assert par.exit_code == 0
    assert sorted(par.stderr.splitlines()) == sorted(seq.stderr.splitlines())


def test_threads_4_stdout_identical_to_threads_1(datadir: Path) -> None:
    """Given six inputs, When screened --threads 1 vs --threads 4 (jobs=1),
    Then stdout is byte-identical (blastn -num_threads never changes results)."""
    files = [str(path) for path in SIX_FILES]
    single = screen(datadir, "--threads", "1", "--nopath", *files)
    multi = screen(datadir, "--threads", "4", "--nopath", *files)
    assert single.exit_code == 0
    assert multi.exit_code == 0
    assert multi.stdout == single.stdout


def test_invalid_jobs_exits_2(datadir: Path) -> None:
    """Given --jobs 0, When screening, Then exit 2 with a USAGE_ERROR envelope."""
    result = screen(datadir, "--jobs", "0", str(CONTIGS / "full.fa"))
    assert result.exit_code == 2
    assert json.loads(result.stderr)["code"] == "USAGE_ERROR"


def test_jobs_4_error_envelope_matches_sequential(datadir: Path, tmp_path: Path) -> None:
    """Given a junk file mid-list, When screened --jobs 1 vs --jobs 4, Then
    both fail identically (same exit code, same JSON envelope)."""
    junk = tmp_path / "junk.txt"
    junk.write_text("this is not sequence data at all\n", encoding="utf-8")
    files = [str(CONTIGS / "full.fa"), str(junk), str(CONTIGS / "none.fa")]
    seq = screen(datadir, "--jobs", "1", *files)
    par = screen(datadir, "--jobs", "4", *files)
    assert seq.exit_code == 5
    assert par.exit_code == 5
    # The envelope is stderr's last line (Processing chatter precedes it).
    assert json.loads(par.stderr.splitlines()[-1]) == json.loads(seq.stderr.splitlines()[-1])
