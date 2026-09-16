"""Unit tests for datadir resolution precedence (gapit.config)."""

from pathlib import Path

import pytest

from gapit.config import resolve_datadir
from gapit.errors import DatabaseError


def make_datadir(path: Path) -> Path:
    path.mkdir(parents=True)
    return path


def test_cli_value_wins_over_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cli_dir = make_datadir(tmp_path / "cli")
    env_dir = make_datadir(tmp_path / "env")
    monkeypatch.setenv("GAPIT_DATADIR", str(env_dir))
    assert resolve_datadir(cli_dir) == cli_dir.resolve()


def test_env_var_used_when_no_cli_value(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    env_dir = make_datadir(tmp_path / "env")
    monkeypatch.setenv("GAPIT_DATADIR", str(env_dir))
    assert resolve_datadir(None) == env_dir.resolve()


def test_default_datadir_when_unset(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("GAPIT_DATADIR", raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))
    default_dir = make_datadir(tmp_path / ".local" / "share" / "gapit" / "db")
    assert resolve_datadir(None) == default_dir.resolve()


def test_tilde_in_env_var_is_expanded(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    data_dir = make_datadir(tmp_path / "data")
    monkeypatch.setenv("GAPIT_DATADIR", "~/data")
    assert resolve_datadir(None) == data_dir.resolve()


def test_missing_datadir_raises_database_error(tmp_path: Path) -> None:
    with pytest.raises(DatabaseError) as excinfo:
        resolve_datadir(tmp_path / "nope")
    assert excinfo.value.exit_code == 4
    assert excinfo.value.code == "DATADIR_NOT_FOUND"
