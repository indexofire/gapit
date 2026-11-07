"""Tests for the `gspec schema` introspection command."""

import json

from click.testing import CliRunner

from gspec.cli import app

runner = CliRunner()


def test_schema_report() -> None:
    result = runner.invoke(app, ["schema", "report"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert "properties" in payload
    assert "schema" in payload["properties"]
    assert "HitDocument" in payload.get("$defs", {})


def test_schema_list_removed_exits_2() -> None:
    """Given the retired gspec.list/1 document, When introspected, Then the
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
    """Given the floors selector, When introspected, Then the gspec.floors/1
    contract surfaces: schema const, [0, 100]-bounded per-gene values."""
    result = runner.invoke(app, ["schema", "floors"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["properties"]["schema"]["const"] == "gspec.floors/1"
    assert payload["properties"]["genes"]["additionalProperties"]["maximum"] == 100.0
    assert payload["properties"]["genes"]["additionalProperties"]["minimum"] == 0.0


def test_schema_unknown_name_exits_2_with_usage_envelope() -> None:
    result = runner.invoke(app, ["schema", "nope"])
    assert result.exit_code == 2
    last_line = [line for line in result.stderr.splitlines() if line.strip()][-1]
    payload = json.loads(last_line)
    assert payload["code"] == "USAGE_ERROR"


def test_schema_cluster_omits_internal_winner_field() -> None:
    """The in-memory PhenotypeDetail.winner is never serialized, so it must
    not be advertised by the public cluster schema."""
    result = runner.invoke(app, ["schema", "cluster"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    detail = payload["$defs"]["PhenotypeDetail"]
    assert "winner" not in detail["properties"]


def test_schema_registers_database_listing_documents() -> None:
    """The database listing documents are first-class agent contracts now."""
    for name, schema in (
        ("dblist", "gspec.dblist/1"),
        ("dbsetup", "gspec.dbsetup/1"),
        ("dboutdated", "gspec.dboutdated/1"),
        ("datadir", "gspec.datadir/1"),
    ):
        result = runner.invoke(app, ["schema", name])
        assert result.exit_code == 0, name
        payload = json.loads(result.stdout)
        assert payload["properties"]["schema"]["const"] == schema


def test_schema_locus_map_is_registered() -> None:
    """Given the rich locus-map document, When introspected, Then its schema
    is a first-class gspec.locus_map/1 contract."""
    result = runner.invoke(app, ["schema", "locus_map"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["properties"]["schema"]["const"] == "gspec.locus_map/1"
    assert "LocusMap" in payload.get("$defs", {})


def test_schema_serotype_is_registered() -> None:
    """Given the serotype profile document, When introspected, Then its
    schema is the gspec.serotype/1 contract."""
    result = runner.invoke(app, ["schema", "serotype"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["properties"]["schema"]["const"] == "gspec.serotype/1"
    assert "SerotypeRule" in payload.get("$defs", {})
