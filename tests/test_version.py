"""Phase 0 smoke tests: package version and CLI entrypoint."""

from typer.testing import CliRunner

from gapit import __version__
from gapit.cli import app

runner = CliRunner()


def test_version_is_semver() -> None:
    major, minor, patch = __version__.split(".")
    assert major.isdigit() and minor.isdigit() and patch.isdigit()


def test_cli_version_flag() -> None:
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.stdout
