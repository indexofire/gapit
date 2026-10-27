"""Integration tests: gapit setupdb against a real fixture datadir.

Exercises the real makeblastdb/blastdbcmd binaries from the pixi environment;
each test copies the committed fixture into its own tmp datadir (isolated).
"""

import re
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from gapit.cli import app
from gapit.db import DatabaseInfo, list_databases

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"

runner = CliRunner()


def make_datadir(tmp_path: Path) -> Path:
    datadir = tmp_path / "datadir"
    shutil.copytree(FIXTURE_DB_DIR, datadir)
    return datadir


def test_setupdb_then_listing_roundtrip(tmp_path: Path) -> None:
    """Given the fixture, When setupdb runs and then databases are listed,
    Then exit 0, all eight bundled databases materialize (one stderr note
    each, alphabetical — vp_ok's counting loci) alongside tinyamr's Indexed
    line, and the listing carries all nine databases."""
    datadir = make_datadir(tmp_path)
    setup = runner.invoke(app, ["setupdb", "--datadir", str(datadir)])
    assert setup.exit_code == 0
    stderr_lines = setup.stderr.splitlines()
    assert stderr_lines[:6] == [
        f"gapit: materializing bundled database {name} ({records} records) into {datadir}"
        for name, records in (
            ("cje_cps", 37),
            ("ecoh", 597),
            ("ecoli_dec", 16),
            ("lm_doumith", 5),
            ("ncbi", 8373),
            ("resfinder", 3206),
        )
    ]
    assert stderr_lines[6:8] == [
        f"gapit: materializing bundled database vc_serogroup (3 records) into {datadir}",
        f"gapit: materializing bundled database vp_ok (195 loci) into {datadir}",
    ]
    assert stderr_lines[8:] == [
        "Indexed cje_cps (37 sequences, nucl)",
        "Indexed ecoh (597 sequences, nucl)",
        "Indexed ecoli_dec (16 sequences, nucl)",
        "Indexed lm_doumith (5 sequences, nucl)",
        "Indexed ncbi (8373 sequences, nucl)",
        "Indexed resfinder (3206 sequences, nucl)",
        "Indexed tinyamr (3 sequences, nucl)",
        "Indexed vc_serogroup (3 sequences, nucl)",
        "Indexed vp_ok (195 sequences, nucl)",
    ]

    infos = list_databases(datadir, setupdb=False)
    assert [info.name for info in infos] == [
        "cje_cps",
        "ecoh",
        "ecoli_dec",
        "lm_doumith",
        "ncbi",
        "resfinder",
        "tinyamr",
        "vc_serogroup",
        "vp_ok",
    ]
    tinyamr = infos[6]
    assert tinyamr.n_sequences == 3
    assert tinyamr.dbtype == "nucl"
    assert re.fullmatch(r"\d{4}-[A-Za-z]{3}-\d{2}", tinyamr.date)


def test_setupdb_on_empty_datadir_materializes_bundled(tmp_path: Path) -> None:
    """Given an empty datadir, When setupdb runs, Then exit 0 with the
    bundled ecoli_dec materialized into it (manifest, typing, BLAST index)."""
    datadir = tmp_path / "empty"
    datadir.mkdir()
    setup = runner.invoke(app, ["setupdb", "--datadir", str(datadir)])
    assert setup.exit_code == 0
    assert "materializing bundled database ecoli_dec" in setup.stderr
    assert "Indexed ecoli_dec (16 sequences, nucl)" in setup.stderr
    assert (datadir / "ecoli_dec" / "gapit-manifest.json").is_file()
    assert (datadir / "ecoli_dec" / "typing.json").is_file()


def test_setupdb_bootstraps_a_missing_datadir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a datadir path that does not exist, When setupdb runs, Then the
    root is created before the bundled/index work (fresh-machine bootstrap)."""
    import gapit.cli as cli_module

    # Keep this a bootstrap test: no bundle materialization or index build.
    def no_bundles() -> list[str]:
        return []

    def no_databases(datadir: Path, *, setupdb: bool, debug: bool = False) -> list[DatabaseInfo]:
        return []

    monkeypatch.setattr(cli_module, "bundled_names", no_bundles)
    monkeypatch.setattr(cli_module.db, "list_databases", no_databases)
    missing = tmp_path / "fresh" / "machine" / "db"

    result = runner.invoke(app, ["setupdb", "--datadir", str(missing)])

    assert result.exit_code == 0, result.stderr
    assert missing.is_dir()


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
    assert "Indexed agtcprot (2 sequences, prot)" in setup.stderr.splitlines()
    assert (db_dir / "sequences.pin").is_file()
    assert not (db_dir / "sequences.nin").exists()
    agtcprot = next(
        info for info in list_databases(datadir, setupdb=False) if info.name == "agtcprot"
    )
    assert agtcprot.dbtype == "prot"
