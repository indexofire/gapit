"""CLI tests for `gapit db path`: resolved datadir introspection."""

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from gapit.cli import app

runner = CliRunner()


def test_db_path_prints_the_resolved_path_without_creating_it(tmp_path: Path) -> None:
    """Given an explicit --datadir, When `db path` runs, Then it prints the
    resolved absolute path and leaves the directory untouched (introspection,
    not bootstrap)."""
    missing = tmp_path / "nested" / "db"

    result = runner.invoke(app, ["db", "path", "--datadir", str(missing)])

    assert result.exit_code == 0, result.stderr
    assert result.stdout.strip() == str(missing.resolve())
    assert not missing.exists()


def test_db_path_json_reports_source_and_existence(tmp_path: Path) -> None:
    """Given a missing explicit datadir, When `db path --json` runs, Then the
    gapit.datadir/1 document names the path, source=cli, and exists=false."""
    missing = tmp_path / "nested" / "db"

    result = runner.invoke(app, ["db", "path", "--datadir", str(missing), "--json"])

    assert result.exit_code == 0, result.stderr
    document = json.loads(result.stdout)
    assert document["schema"] == "gapit.datadir/1"
    assert document["datadir"] == str(missing.resolve())
    assert document["source"] == "cli"
    assert document["exists"] is False
    assert not missing.exists()


def test_db_path_uses_xdg_data_home(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given XDG_DATA_HOME and no GAPIT_DATADIR/--datadir, When `db path`
    runs, Then the XDG user-data root is selected and reported as source=xdg."""
    monkeypatch.delenv("GAPIT_DATADIR", raising=False)
    xdg = tmp_path / "xdg"
    monkeypatch.setenv("XDG_DATA_HOME", str(xdg))

    result = runner.invoke(app, ["db", "path", "--json"])

    assert result.exit_code == 0, result.stderr
    document = json.loads(result.stdout)
    assert document["datadir"] == str((xdg / "gapit" / "db").resolve())
    assert document["source"] == "xdg"
    assert document["exists"] is False


def test_db_path_prefers_gapit_datadir_over_xdg(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GAPIT_DATADIR", str(tmp_path / "env"))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))

    result = runner.invoke(app, ["db", "path", "--json"])

    assert result.exit_code == 0, result.stderr
    document = json.loads(result.stdout)
    assert document["datadir"] == str((tmp_path / "env").resolve())
    assert document["source"] == "env"
