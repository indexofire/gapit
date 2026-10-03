"""Tests for the `gapit schema` introspection command."""

import json

from typer.testing import CliRunner

from gapit.cli import app

runner = CliRunner()


def test_schema_report() -> None:
    result = runner.invoke(app, ["schema", "report"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert "properties" in payload
    assert "schema" in payload["properties"]
    assert "HitDocument" in payload.get("$defs", {})


def test_schema_list_removed_exits_2() -> None:
    """Given the retired gapit.list/1 document, When introspected, Then the
    name is rejected as unknown with a usage envelope on stderr."""
    result = runner.invoke(app, ["schema", "list"])
    assert result.exit_code == 2
    assert "unknown schema name: list" in result.stderr


def test_schema_error() -> None:
    result = runner.invoke(app, ["schema", "error"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    for key in ("schema", "code", "message", "context"):
        assert key in payload["properties"]


def test_schema_version() -> None:
    result = runner.invoke(app, ["schema", "version"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert "version" in payload["properties"]


def test_schema_floors() -> None:
    """Given the floors selector, When introspected, Then the gapit.floors/1
    contract surfaces: schema const, [0, 100]-bounded per-gene values."""
    result = runner.invoke(app, ["schema", "floors"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["properties"]["schema"]["const"] == "gapit.floors/1"
    assert payload["properties"]["genes"]["additionalProperties"]["maximum"] == 100.0
    assert payload["properties"]["genes"]["additionalProperties"]["minimum"] == 0.0


def test_schema_unknown_name_exits_2_with_usage_envelope() -> None:
    result = runner.invoke(app, ["schema", "nope"])
    assert result.exit_code == 2
    last_line = [line for line in result.stderr.splitlines() if line.strip()][-1]
    payload = json.loads(last_line)
    assert payload["code"] == "USAGE_ERROR"
