"""Phase 0 smoke tests: package version and CLI entrypoint."""

import json
import tomllib
from pathlib import Path

from click.testing import CliRunner

from gapit import __version__
from gapit.cli import app

runner = CliRunner()


def test_version_is_semver() -> None:
    major, minor, patch = __version__.split(".")
    assert major.isdigit() and minor.isdigit() and patch.isdigit()


def test_version_matches_pyproject() -> None:
    """Given the source tree, When read, Then __version__ matches the
    version pyproject.toml ships — release bumps cannot silently drift."""
    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    document = tomllib.loads(pyproject.read_text(encoding="utf-8"))
    assert document["project"]["version"] == __version__


def test_cli_version_flag() -> None:
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert result.stdout == f"gapit {__version__}\n"


def test_cli_version_json() -> None:
    """Given --version --json, When run, Then one compact line with the
    gapit.version/1 contract."""
    result = runner.invoke(app, ["--version", "--json"])
    assert result.exit_code == 0
    assert "\n" not in result.stdout.rstrip("\n")
    payload = json.loads(result.stdout)
    assert payload == {
        "schema": "gapit.version/1",
        "name": "gapit",
        "version": __version__,
    }
