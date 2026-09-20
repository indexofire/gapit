"""Pipeline tests for blast's native-normalization stage: the normalized FASTA
must actually reach blastn on stdin, blast failures carry their stderr as
BLAST_FAILED, and a missing blastn is a typed dependency error. (The retired
``any2fasta | blastn`` Popen chain needed a stderr-drain thread; the
``subprocess.run(input=...)`` rewrite makes that deadlock structurally
impossible, so fakes here only stand in for blastn.) Fakes are executable
python scripts named ``blastn`` on a monkeypatched PATH.
"""

import os
from pathlib import Path

import pytest

from gapit.blast import run_screen
from gapit.db import Database
from gapit.errors import DependencyError, GapitError
from gapit.report import ScreeningParams

LOWERCASE_FASTA = ">contig1\nacgtacgtac\n"


@pytest.fixture()
def fakedb(tmp_path: Path) -> Database:
    """A Database whose paths are never read (the fake binaries ignore argv)."""
    return Database(name="fakedb", path=tmp_path, sequences_path=tmp_path / "sequences")


def fake_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")
    return bin_dir


def write_fake_blastn(bin_dir: Path, script_body: str) -> None:
    binary = bin_dir / "blastn"
    binary.write_text("#!/usr/bin/env python3\n" + script_body, encoding="utf-8")
    binary.chmod(0o755)


def test_blastn_receives_normalized_fasta_on_stdin(
    tmp_path: Path, fakedb: Database, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a lowercase FASTA query, When screened through a fake blastn that
    dumps stdin to a file, Then the dump is the uppercase re-emitted FASTA —
    the native converter feeds blast verbatim (the old pipe's contract)."""
    bin_dir = fake_path(monkeypatch, tmp_path)
    dump = tmp_path / "stdin.txt"
    monkeypatch.setenv("GAPIT_STDIN_DUMP", str(dump))
    write_fake_blastn(
        bin_dir,
        "import os, sys\n"
        "with open(os.environ['GAPIT_STDIN_DUMP'], 'w') as fh:\n"
        "    fh.write(sys.stdin.read())\n"
        "sys.exit(0)\n",
    )
    query = tmp_path / "contigs.fa"
    query.write_text(LOWERCASE_FASTA, encoding="utf-8")
    params = ScreeningParams(db="fakedb", minid=80.0, mincov=80.0, threads=1)
    report = run_screen(query, fakedb, params, dbtype="nucl")
    assert report == []
    assert dump.read_text(encoding="utf-8") == ">contig1\nACGTACGTAC\n"


def test_blastn_failure_raises_blast_failed_with_stderr(
    tmp_path: Path, fakedb: Database, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a blastn that exits 1 with a stderr message, When screened, Then
    GapitError BLAST_FAILED carries that stderr text."""
    bin_dir = fake_path(monkeypatch, tmp_path)
    write_fake_blastn(
        bin_dir,
        "import sys\nsys.stderr.write('kaboom: fake blast failure\\n')\nsys.exit(1)\n",
    )
    query = tmp_path / "contigs.fa"
    query.write_text(LOWERCASE_FASTA, encoding="utf-8")
    params = ScreeningParams(db="fakedb", minid=80.0, mincov=80.0, threads=1)
    with pytest.raises(GapitError) as excinfo:
        run_screen(query, fakedb, params, dbtype="nucl")
    assert excinfo.value.code == "BLAST_FAILED"
    assert "kaboom" in str(excinfo.value)
    assert excinfo.value.context["binary"] == "blastn"


def test_missing_blastn_raises_dependency_error(
    tmp_path: Path, fakedb: Database, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a PATH without blastn, When screened, Then DependencyError
    MISSING_DEPENDENCY names blastn (exit 3)."""
    empty_bin = tmp_path / "empty-bin"
    empty_bin.mkdir()
    monkeypatch.setenv("PATH", str(empty_bin))
    query = tmp_path / "contigs.fa"
    query.write_text(LOWERCASE_FASTA, encoding="utf-8")
    params = ScreeningParams(db="fakedb", minid=80.0, mincov=80.0, threads=1)
    with pytest.raises(DependencyError) as excinfo:
        run_screen(query, fakedb, params, dbtype="nucl")
    assert excinfo.value.code == "MISSING_DEPENDENCY"
    assert excinfo.value.context["binary"] == "blastn"
