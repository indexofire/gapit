"""Unit tests for datadir resolution precedence (gapit.config)."""

import stat
from pathlib import Path

import pytest

from gapit.config import datadir_path, datadir_source, ensure_datadir, resolve_datadir
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
    monkeypatch.delenv("XDG_DATA_HOME", raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))
    default_dir = make_datadir(tmp_path / ".local" / "share" / "gapit" / "db")
    assert resolve_datadir(None) == default_dir.resolve()


def test_xdg_data_home_used_when_unset(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given XDG_DATA_HOME and no GAPIT_DATADIR, When resolved, Then the XDG
    user-data root wins over the hardcoded ~/.local/share fallback."""
    monkeypatch.delenv("GAPIT_DATADIR", raising=False)
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    path, source = datadir_source(None)
    assert path == (tmp_path / "xdg" / "gapit" / "db").resolve()
    assert source == "xdg"


def test_gapit_datadir_wins_over_xdg(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GAPIT_DATADIR", str(tmp_path / "env"))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    path, source = datadir_source(None)
    assert path == (tmp_path / "env").resolve()
    assert source == "env"


def test_cli_datadir_wins_and_is_reported(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GAPIT_DATADIR", str(tmp_path / "env"))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    path, source = datadir_source(tmp_path / "cli")
    assert path == (tmp_path / "cli").resolve()
    assert source == "cli"


def test_relative_xdg_data_home_is_ignored(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The XDG spec requires an absolute XDG_DATA_HOME; a relative value is
    ignored and the default fallback is used."""
    monkeypatch.delenv("GAPIT_DATADIR", raising=False)
    monkeypatch.setenv("XDG_DATA_HOME", "relative/data")
    monkeypatch.setenv("HOME", str(tmp_path))
    path, source = datadir_source(None)
    assert path == (tmp_path / ".local" / "share" / "gapit" / "db").resolve()
    assert source == "default"


def test_empty_gapit_datadir_falls_back_to_xdg(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("GAPIT_DATADIR", "")
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    path, source = datadir_source(None)
    assert path == (tmp_path / "xdg" / "gapit" / "db").resolve()
    assert source == "xdg"


def test_datadir_path_does_not_require_existence(tmp_path: Path) -> None:
    missing = tmp_path / "missing" / "db"
    assert datadir_path(missing) == missing.resolve()
    assert not missing.exists()


def test_ensure_datadir_creates_a_private_root(tmp_path: Path) -> None:
    missing = tmp_path / "private" / "db"
    assert ensure_datadir(missing) == missing.resolve()
    assert missing.is_dir()
    assert stat.S_IMODE(missing.stat().st_mode) == 0o700


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
