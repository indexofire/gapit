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


def test_schema_list() -> None:
    result = runner.invoke(app, ["schema", "list"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert "databases" in payload["properties"]


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


def test_schema_unknown_name_exits_2_with_usage_envelope() -> None:
    result = runner.invoke(app, ["schema", "nope"])
    assert result.exit_code == 2
    last_line = [line for line in result.stderr.splitlines() if line.strip()][-1]
    payload = json.loads(last_line)
    assert payload["code"] == "USAGE_ERROR"
