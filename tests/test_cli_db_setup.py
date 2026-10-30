"""Integration tests: gapit db setup against a real fixture datadir.

Exercises the real makeblastdb/blastdbcmd binaries from the pixi environment;
each test copies the committed fixture into its own tmp datadir (isolated).
"""

import json
import re
import shutil
from collections.abc import Iterable
from pathlib import Path

import pytest
from click.testing import CliRunner

from gapit.cli import app
from gapit.db import DatabaseInfo, list_databases
from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"

runner = CliRunner()


def make_datadir(tmp_path: Path) -> Path:
    datadir = tmp_path / "datadir"
    shutil.copytree(FIXTURE_DB_DIR, datadir)
    return datadir


def patch_synthetic_provider(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, str]:
    """Point both registry views at one file:// provider and create a datadir."""
    name = "synth"
    source = tmp_path / "upstream.fa"
    source.write_text(">syn_a first synthetic gene\nACGTAGACGTAG\n", encoding="utf-8")

    def transform(workdir: Path) -> Iterable[Record]:
        for record in iter_fasta(workdir / "upstream.fa"):
            yield Record(db=name, gene=record.id, sequence=record.sequence)

    provider = Provider(
        name=name,
        description="synthetic provider for db setup",
        vendor="Synthetica",
        source_urls=(source.as_uri(),),
        dbtype="nucl",
        transform=transform,
    )
    registry = {name: provider}
    # perform_setup reads gapit.providers.REGISTRY for membership and
    # perform_fetch reads gapit.db_ops.REGISTRY for the download; both must
    # see the synthetic catalogue.
    monkeypatch.setattr("gapit.providers.REGISTRY", registry)
    monkeypatch.setattr("gapit.db_ops.REGISTRY", registry)
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    return datadir, name


def test_setup_then_listing_roundtrip(tmp_path: Path) -> None:
    """Given the fixture, When db setup runs and then databases are listed,
    Then exit 0, all four bundled databases materialize (one stderr note
    each, alphabetical — vpa_ok's counting loci) alongside tinyamr's Indexed
    line, and the listing carries all five databases."""
    datadir = make_datadir(tmp_path)
    setup = runner.invoke(app, ["db", "setup", "--datadir", str(datadir)])
    assert setup.exit_code == 0
    stderr_lines = setup.stderr.splitlines()
    assert stderr_lines[:3] == [
        f"gapit: materializing bundled database {name} ({records} records) into {datadir}"
        for name, records in (
            ("cje_cps", 42),
            ("ecoli_dec", 16),
            ("lmo_sero", 5),
        )
    ]
    assert stderr_lines[3:4] == [
        f"gapit: materializing bundled database vpa_ok (195 loci) into {datadir}",
    ]
    assert stderr_lines[4:] == [
        "gapit: indexed cje_cps (42 sequences, nucl)",
        "gapit: indexed ecoli_dec (16 sequences, nucl)",
        "gapit: indexed lmo_sero (5 sequences, nucl)",
        "gapit: indexed tinyamr (3 sequences, nucl)",
        "gapit: indexed vpa_ok (195 sequences, nucl)",
    ]

    infos = list_databases(datadir, reindex=False)
    assert [info.name for info in infos] == [
        "cje_cps",
        "ecoli_dec",
        "lmo_sero",
        "tinyamr",
        "vpa_ok",
    ]
    tinyamr = infos[3]
    assert tinyamr.n_sequences == 3
    assert tinyamr.dbtype == "nucl"
    assert re.fullmatch(r"\d{4}-[A-Za-z]{3}-\d{2}", tinyamr.date)


def test_setup_on_empty_datadir_materializes_bundled(tmp_path: Path) -> None:
    """Given an empty datadir, When db setup runs, Then exit 0 with the
    bundled ecoli_dec materialized into it (manifest, typing, BLAST index)."""
    datadir = tmp_path / "empty"
    datadir.mkdir()
    setup = runner.invoke(app, ["db", "setup", "--datadir", str(datadir)])
    assert setup.exit_code == 0
    assert "materializing bundled database ecoli_dec" in setup.stderr
    assert "gapit: indexed ecoli_dec (16 sequences, nucl)" in setup.stderr
    assert (datadir / "ecoli_dec" / "gapit-manifest.json").is_file()
    assert (datadir / "ecoli_dec" / "typing.json").is_file()


def test_setup_bootstraps_a_missing_datadir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a datadir path that does not exist, When db setup runs, Then the
    root is created before the bundled/index work (fresh-machine bootstrap)."""
    import gapit.db_setup_ops as setup_module

    # Keep this a bootstrap test: no bundle materialization or index build.
    def no_bundles() -> list[str]:
        return []

    def no_databases(datadir: Path, *, reindex: bool, debug: bool = False) -> list[DatabaseInfo]:
        return []

    monkeypatch.setattr(setup_module, "bundled_names", no_bundles)
    monkeypatch.setattr(setup_module.db, "list_databases", no_databases)
    missing = tmp_path / "fresh" / "machine" / "db"

    result = runner.invoke(app, ["db", "setup", "--datadir", str(missing)])

    assert result.exit_code == 0, result.stderr
    assert missing.is_dir()


def test_setup_debug_echoes_makeblastdb_argv(tmp_path: Path) -> None:
    """Given an unindexed fixture datadir, When `db setup --debug` runs, Then
    exit 0 and the makeblastdb argv is echoed to stderr as a `gapit: run:`
    line (the no-flag default keeps stderr to the Indexed line only — pinned
    by test_setup_then_listing_roundtrip)."""
    datadir = make_datadir(tmp_path)
    result = runner.invoke(app, ["db", "setup", "--debug", "--datadir", str(datadir)])
    assert result.exit_code == 0
    run_lines = [line for line in result.stderr.splitlines() if line.startswith("gapit: run:")]
    assert run_lines and run_lines[0].startswith("gapit: run: makeblastdb -in ")


def test_setup_honors_manifest_dbtype_on_reindex(tmp_path: Path) -> None:
    """Given a gapit-built prot database whose sequences are pure A/G/T/C
    (the abricate mol_type heuristic alone would say nucl), When the BLAST
    index is deleted and rebuilt via db setup, Then the manifest dbtype wins:
    a .pin index exists, .nin does not, and list_databases reports prot.
    Manifest-less (abricate-built) dirs keep the heuristic — pinned by
    test_setup_then_listing_roundtrip."""
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
    setup = runner.invoke(app, ["db", "setup", "--datadir", str(datadir)])
    assert setup.exit_code == 0
    assert "gapit: indexed agtcprot (2 sequences, prot)" in setup.stderr.splitlines()
    assert (db_dir / "sequences.pin").is_file()
    assert not (db_dir / "sequences.nin").exists()
    agtcprot = next(
        info for info in list_databases(datadir, reindex=False) if info.name == "agtcprot"
    )
    assert agtcprot.dbtype == "prot"


def test_setup_json_receipt_on_existing_datadir(tmp_path: Path) -> None:
    """Given an installed datadir, When `db setup --json --no-bundles` runs,
    Then stdout is one gapit.dbsetup/1 document (no progress noise) listing
    the indexed databases and an empty materialized list."""
    datadir = make_datadir(tmp_path)
    result = runner.invoke(
        app, ["db", "setup", "--json", "--no-bundles", "--datadir", str(datadir)]
    )
    assert result.exit_code == 0
    assert result.stderr == ""
    document = json.loads(result.stdout)
    assert document["schema"] == "gapit.dbsetup/1"
    assert document["datadir"] == str(datadir)
    assert document["materialized"] == []
    assert [entry["db"] for entry in document["databases"]] == ["tinyamr"]
    assert document["databases"][0]["dbtype"] == "nucl"


def test_setup_no_bundles_on_empty_datadir_leaves_it_empty(tmp_path: Path) -> None:
    """Given an empty datadir, When `db setup --no-bundles` runs, Then the
    root is created but no bundled database is materialized."""
    datadir = tmp_path / "empty"
    datadir.mkdir()
    result = runner.invoke(
        app, ["db", "setup", "--json", "--no-bundles", "--datadir", str(datadir)]
    )
    assert result.exit_code == 0
    document = json.loads(result.stdout)
    assert document["materialized"] == []
    assert document["databases"] == []
    assert sorted(path.name for path in datadir.iterdir()) == []


def test_setup_quiet_suppresses_progress(tmp_path: Path) -> None:
    """Given --quiet, When db setup runs, Then stderr carries no progress."""
    datadir = make_datadir(tmp_path)
    result = runner.invoke(app, ["db", "setup", "-q", "--datadir", str(datadir)])
    assert result.exit_code == 0
    assert result.stderr == ""


def test_old_setupdb_spelling_is_removed() -> None:
    """Given either former spelling, When invoked, Then click reports an
    unknown command: the alias is gone, not merely hidden."""
    for argv in (["setupdb"], ["db", "setupdb"]):
        result = runner.invoke(app, argv)
        assert result.exit_code == 2, argv
        assert "No such command 'setupdb'" in result.stderr, argv


def test_setup_named_database_indexes_only_that(tmp_path: Path) -> None:
    """Given an installed datadir, When db setup names one database, Then
    only that database is indexed and no bundle is materialized."""
    datadir = make_datadir(tmp_path)
    result = runner.invoke(
        app, ["db", "setup", "-d", "tinyamr", "--json", "--datadir", str(datadir)]
    )
    assert result.exit_code == 0, result.stderr
    document = json.loads(result.stdout)
    assert document["materialized"] == []
    assert [entry["db"] for entry in document["databases"]] == ["tinyamr"]
    assert not (datadir / "cje_cps").exists()


def test_setup_named_bundled_database_materializes_only_that(tmp_path: Path) -> None:
    """Given an empty datadir, When db setup names one bundled database,
    Then only that bundle is materialized and indexed."""
    datadir = tmp_path / "empty"
    datadir.mkdir()
    result = runner.invoke(
        app, ["db", "setup", "-d", "ecoli_dec", "--json", "--datadir", str(datadir)]
    )
    assert result.exit_code == 0, result.stderr
    document = json.loads(result.stdout)
    assert document["materialized"] == ["ecoli_dec"]
    assert [entry["db"] for entry in document["databases"]] == ["ecoli_dec"]
    visible = sorted(
        path.name for path in datadir.iterdir() if not path.name.startswith(".gapit-txn")
    )
    assert visible == ["ecoli_dec"]


def test_setup_unknown_name_is_database_not_found(tmp_path: Path) -> None:
    """Given a name that is neither installed nor bundled, When db setup runs,
    Then it exits 4 with DATABASE_NOT_FOUND before writing anything."""
    datadir = tmp_path / "empty"
    datadir.mkdir()
    result = runner.invoke(
        app, ["db", "setup", "-d", "nosuch", "--json", "--datadir", str(datadir)]
    )
    assert result.exit_code == 4
    envelope = json.loads(result.stderr.splitlines()[-1])
    assert envelope["code"] == "DATABASE_NOT_FOUND"
    assert "nosuch" in envelope["message"]
    assert sorted(path.name for path in datadir.iterdir()) == []


def test_setup_no_bundles_refuses_a_missing_bundled_name(tmp_path: Path) -> None:
    """Given --no-bundles and an uninstalled bundled name, When db setup runs,
    Then it exits 4 instead of silently doing nothing."""
    datadir = tmp_path / "empty"
    datadir.mkdir()
    result = runner.invoke(app, ["db", "setup", "-d", "ecoli_dec", "-B", "--datadir", str(datadir)])
    assert result.exit_code == 4
    envelope = json.loads(result.stderr.splitlines()[-1])
    assert envelope["code"] == "DATABASE_NOT_FOUND"
    assert "--no-bundles" in envelope["message"]


def test_setup_all_flag_indexes_installed_without_bundles(tmp_path: Path) -> None:
    """Given an installed datadir, When db setup uses --all --no-bundles, Then
    it indexes the installed database and materializes nothing."""
    datadir = make_datadir(tmp_path)
    result = runner.invoke(
        app,
        ["db", "setup", "-a", "-B", "--json", "--datadir", str(datadir)],
    )
    assert result.exit_code == 0, result.stderr
    document = json.loads(result.stdout)
    assert document["materialized"] == []
    assert document["fetched"] == []
    assert [entry["db"] for entry in document["databases"]] == ["tinyamr"]


def test_setup_all_and_db_are_mutually_exclusive(tmp_path: Path) -> None:
    """Given --all and --db together, When db setup runs, Then it is a usage
    error before any datadir work."""
    result = runner.invoke(
        app,
        ["db", "setup", "-a", "-d", "tinyamr", "--datadir", str(tmp_path / "db")],
    )
    assert result.exit_code == 2
    assert "not both" in result.stderr


def test_setup_named_provider_fetches_then_indexes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a named provider that is not installed, When db setup -d runs,
    Then it fetches the provider and indexes it in the same call."""
    datadir, name = patch_synthetic_provider(tmp_path, monkeypatch)
    result = runner.invoke(app, ["db", "setup", "-d", name, "--json", "--datadir", str(datadir)])
    assert result.exit_code == 0, result.stderr
    document = json.loads(result.stdout)
    assert document["materialized"] == []
    assert document["fetched"] == [name]
    assert [entry["db"] for entry in document["databases"]] == [name]
    assert (datadir / name / "sequences").is_file()
