"""Integration tests: gapit list / setupdb against a real fixture datadir.

Exercises the real makeblastdb/blastdbcmd binaries from the pixi environment;
each test copies the committed fixture into its own tmp datadir (isolated).
"""

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
    """Given an unindexed fixture datadir, When listed, Then exit 4 with an
    ERROR line on stderr and no table on stdout."""
    datadir = make_datadir(tmp_path)
    result = runner.invoke(app, ["list", "--datadir", str(datadir)])
    assert result.exit_code == 4
    assert result.stdout == ""
    assert "ERROR:" in result.stderr
    assert "not indexed" in result.stderr


def test_list_with_missing_datadir_exits_4(tmp_path: Path) -> None:
    result = runner.invoke(app, ["list", "--datadir", str(tmp_path / "nope")])
    assert result.exit_code == 4
    assert "ERROR:" in result.stderr


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
