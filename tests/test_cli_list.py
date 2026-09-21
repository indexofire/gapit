"""Integration tests: gapit list / setupdb against a real fixture datadir.

Exercises the real makeblastdb/blastdbcmd binaries from the pixi environment;
each test copies the committed fixture into its own tmp datadir (isolated).
"""

import json
import re
import shutil
from pathlib import Path

from pydantic import BaseModel, Field, TypeAdapter
from typer.testing import CliRunner

from gapit.cli import app

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"

runner = CliRunner()


class ListEntry(BaseModel):
    """One entry of the gapit.list/1 payload (boundary parse for this test)."""

    name: str
    sequences: int
    dbtype: str
    date: str


class ListPayload(BaseModel):
    """The gapit.list/1 payload (`schema` clashes with BaseModel.schema, hence the alias)."""

    schema_name: str = Field(alias="schema")
    databases: list[ListEntry]


def make_datadir(tmp_path: Path) -> Path:
    datadir = tmp_path / "datadir"
    shutil.copytree(FIXTURE_DB_DIR, datadir)
    return datadir


def test_list_before_setupdb_exits_4(tmp_path: Path) -> None:
    """Given an unindexed fixture datadir, When listed, Then exit 4 with a
    gapit.error/1 envelope on stderr and no table on stdout."""
    datadir = make_datadir(tmp_path)
    result = runner.invoke(app, ["list", "--datadir", str(datadir)])
    assert result.exit_code == 4
    assert result.stdout == ""
    envelope = json.loads(result.stderr)
    assert envelope["code"] == "DATABASE_NOT_INDEXED"
    assert "not indexed" in envelope["message"]


def test_list_with_missing_datadir_exits_4(tmp_path: Path) -> None:
    result = runner.invoke(app, ["list", "--datadir", str(tmp_path / "nope")])
    assert result.exit_code == 4
    envelope = json.loads(result.stderr)
    assert envelope["code"] == "DATADIR_NOT_FOUND"


def test_setupdb_then_list_roundtrip(tmp_path: Path) -> None:
    """Given the fixture, When setupdb runs and then list, Then exit 0, one
    Indexed stderr line, and an abricate-style table row for tinyamr."""
    datadir = make_datadir(tmp_path)
    setup = runner.invoke(app, ["setupdb", "--datadir", str(datadir)])
    assert setup.exit_code == 0
    assert setup.stderr.strip() == "Indexed tinyamr (3 sequences, nucl)"

    result = runner.invoke(app, ["list", "--datadir", str(datadir)])
    assert result.exit_code == 0
    lines = result.stdout.splitlines()
    assert lines[0] == "DATABASE\tSEQUENCES\tDBTYPE\tDATE"
    assert len(lines) == 2
    name, sequences, dbtype, date = lines[1].split("\t")
    assert name == "tinyamr"
    assert sequences == "3"
    assert dbtype == "nucl"
    assert re.fullmatch(r"\d{4}-[A-Za-z]{3}-\d{2}", date)


def test_setupdb_runs_list_without_list(tmp_path: Path) -> None:
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
    by test_setupdb_then_list_roundtrip)."""
    datadir = make_datadir(tmp_path)
    result = runner.invoke(app, ["setupdb", "--debug", "--datadir", str(datadir)])
    assert result.exit_code == 0
    run_lines = [line for line in result.stderr.splitlines() if line.startswith("gapit: run:")]
    assert run_lines and run_lines[0].startswith("gapit: run: makeblastdb -in ")


def test_setupdb_honors_manifest_dbtype_on_reindex(tmp_path: Path) -> None:
    """Given a gapit-built prot database whose sequences are pure A/G/T/C
    (the abricate mol_type heuristic alone would say nucl), When the BLAST
    index is deleted and rebuilt via setupdb, Then the manifest dbtype wins:
    a .pin index exists, .nin does not, and list reports prot. Manifest-less
    (abricate-built) dirs keep the heuristic — pinned by
    test_setupdb_then_list_roundtrip."""
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
    result = runner.invoke(app, ["list", "--datadir", str(datadir)])
    assert result.exit_code == 0
    lines = result.stdout.splitlines()
    assert len(lines) == 2
    assert lines[1].split("\t")[2] == "prot"


def test_list_json_output(tmp_path: Path) -> None:
    """Given an indexed fixture datadir, When listed with --json, Then the payload
    parses, carries schema gapit.list/1, and keys come in documented order."""
    datadir = make_datadir(tmp_path)
    assert runner.invoke(app, ["setupdb", "--datadir", str(datadir)]).exit_code == 0
    result = runner.invoke(app, ["list", "--datadir", str(datadir), "--json"])
    assert result.exit_code == 0
    payload = TypeAdapter(ListPayload).validate_json(result.stdout)
    assert payload.schema_name == "gapit.list/1"
    (entry,) = payload.databases
    assert entry.name == "tinyamr"
    assert entry.sequences == 3
    assert entry.dbtype == "nucl"
    assert re.fullmatch(r"\d{4}-[A-Za-z]{3}-\d{2}", entry.date)
    key_positions = [result.stdout.index(f'"{key}"') for key in ListEntry.model_fields]
    assert key_positions == sorted(key_positions)
