"""Integration tests: gapit screen CLI over the real pipeline on tinyamr fixtures."""

import shutil
from pathlib import Path

import pytest
from typer.testing import CliRunner, Result

from gapit.cli import app
from gapit.db import make_blast_db

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"
CONTIGS = Path(__file__).parent / "data" / "contigs"
GOLDEN = Path(__file__).parent / "golden"
MULTI_FILES = ["full.fa", "gap.fa", "none.fa", "sort.fa"]

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
    """Given the same run with --csv, When screened, Then stdout is
    byte-identical to the committed golden CSV."""
    result = screen(datadir, "--csv", "--nopath", *[str(CONTIGS / name) for name in MULTI_FILES])
    assert result.exit_code == 0
    assert result.stdout == (GOLDEN / "tinyamr_multi_nopath.csv").read_text(encoding="utf-8")


def test_format_csv_and_csv_flag_coexist(datadir: Path) -> None:
    """Given both --format csv and --csv, When screened, Then csv wins with no
    error and output matches the --csv-only run."""
    both = screen(datadir, "--format", "csv", "--csv", "--nopath", str(CONTIGS / "full.fa"))
    only = screen(datadir, "--csv", "--nopath", str(CONTIGS / "full.fa"))
    assert both.exit_code == 0
    assert both.stdout == only.stdout


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
    assert "ERROR" in result.stderr


def test_invalid_mincov_exits_2(datadir: Path) -> None:
    result = screen(datadir, "--mincov", "101", str(CONTIGS / "full.fa"))
    assert result.exit_code == 2


def test_invalid_threads_exits_2(datadir: Path) -> None:
    result = screen(datadir, "--threads", "0", str(CONTIGS / "full.fa"))
    assert result.exit_code == 2


def test_missing_input_file_exits_5(datadir: Path) -> None:
    result = screen(datadir, str(datadir / "nope.fa"))
    assert result.exit_code == 5
    assert "ERROR" in result.stderr


def test_junk_input_exits_5(datadir: Path, tmp_path: Path) -> None:
    junk = tmp_path / "junk.txt"
    junk.write_text("this is not sequence data at all\n", encoding="utf-8")
    result = screen(datadir, str(junk))
    assert result.exit_code == 5


def test_unknown_db_exits_4_and_lists_available(datadir: Path) -> None:
    result = screen(datadir, "--db", "nope", str(CONTIGS / "full.fa"))
    assert result.exit_code == 4
    assert "tinyamr" in result.stderr
    assert "Available" in result.stderr


def test_no_input_files_exits_2(datadir: Path) -> None:
    result = screen(datadir)
    assert result.exit_code == 2
    assert "ERROR" in result.stderr


def test_debug_echoes_blast_argv(datadir: Path) -> None:
    result = screen(datadir, "--debug", str(CONTIGS / "full.fa"))
    assert result.exit_code == 0
    assert "DEBUG:" in result.stderr
    assert "blastn" in result.stderr
    assert "-perc_identity" in result.stderr
