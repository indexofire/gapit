"""Integration tests: gapit setupdb against a real fixture datadir.

Exercises the real makeblastdb/blastdbcmd binaries from the pixi environment;
each test copies the committed fixture into its own tmp datadir (isolated).
"""

import re
import shutil
from pathlib import Path

from typer.testing import CliRunner

from gapit.cli import app
from gapit.db import list_databases

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"

runner = CliRunner()


def make_datadir(tmp_path: Path) -> Path:
    datadir = tmp_path / "datadir"
    shutil.copytree(FIXTURE_DB_DIR, datadir)
    return datadir


def test_setupdb_then_listing_roundtrip(tmp_path: Path) -> None:
    """Given the fixture, When setupdb runs and then databases are listed,
    Then exit 0, one Indexed stderr line, and an abricate-style row for
    tinyamr."""
    datadir = make_datadir(tmp_path)
    setup = runner.invoke(app, ["setupdb", "--datadir", str(datadir)])
    assert setup.exit_code == 0
    assert setup.stderr.strip() == "Indexed tinyamr (3 sequences, nucl)"

    (info,) = list_databases(datadir, setupdb=False)
    assert info.name == "tinyamr"
    assert info.n_sequences == 3
    assert info.dbtype == "nucl"
    assert re.fullmatch(r"\d{4}-[A-Za-z]{3}-\d{2}", info.date)


def test_setupdb_on_empty_datadir_exits_0(tmp_path: Path) -> None:
    """Given an empty datadir, When setupdb runs, Then exit 0 with nothing indexed."""
    datadir = tmp_path / "empty"
    datadir.mkdir()
    setup = runner.invoke(app, ["setupdb", "--datadir", str(datadir)])
    assert setup.exit_code == 0
    assert setup.stderr == ""


def test_setupdb_debug_echoes_makeblastdb_argv(tmp_path: Path) -> None:
    """Given an unindexed fixture datadir, When `setupdb --debug` runs, Then
    exit 0 and the makeblastdb argv is echoed to stderr as a `gapit: run:`
    line (the no-flag default keeps stderr to the Indexed line only — pinned
    by test_setupdb_then_listing_roundtrip)."""
    datadir = make_datadir(tmp_path)
    result = runner.invoke(app, ["setupdb", "--debug", "--datadir", str(datadir)])
    assert result.exit_code == 0
    run_lines = [line for line in result.stderr.splitlines() if line.startswith("gapit: run:")]
    assert run_lines and run_lines[0].startswith("gapit: run: makeblastdb -in ")


def test_setupdb_honors_manifest_dbtype_on_reindex(tmp_path: Path) -> None:
    """Given a gapit-built prot database whose sequences are pure A/G/T/C
    (the abricate mol_type heuristic alone would say nucl), When the BLAST
    index is deleted and rebuilt via setupdb, Then the manifest dbtype wins:
    a .pin index exists, .nin does not, and list_databases reports prot.
    Manifest-less (abricate-built) dirs keep the heuristic — pinned by
    test_setupdb_then_listing_roundtrip."""
    prot_fa = tmp_path / "agtc_prot.fa"
    prot_fa.write_text(
        ">agtc_strep synthetic AGTC-heavy protein\n" + "AGTCAGTC" * 6 + "\n"
        ">agtc_mix synthetic AGTC-heavy protein\n" + "GATCGATC" * 6 + "\n",
        encoding="utf-8",
    )
    datadir = tmp_path / "datadir"
    build = runner.invoke(
        app,
        ["db", "build", "agtcprot", str(prot_fa), "--dbtype", "prot", "--datadir", str(datadir)],
    )
    assert build.exit_code == 0
    db_dir = datadir / "agtcprot"
    for index_file in db_dir.glob("sequences.[np]??"):
        index_file.unlink()
    setup = runner.invoke(app, ["setupdb", "--datadir", str(datadir)])
    assert setup.exit_code == 0
    assert setup.stderr.strip() == "Indexed agtcprot (2 sequences, prot)"
    assert (db_dir / "sequences.pin").is_file()
    assert not (db_dir / "sequences.nin").exists()
    (info,) = list_databases(datadir, setupdb=False)
    assert info.dbtype == "prot"
