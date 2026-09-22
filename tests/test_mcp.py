"""MCP stdio server tests: JSON-RPC handshake, tool listing, and tool calls
(the screen tool drives the real BLAST pipeline on the tinyamr fixture)."""

import json
import shutil
from io import StringIO
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from gapit.cli import app
from gapit.db import make_blast_db
from gapit.mcp import serve

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"
CONTIGS = Path(__file__).parent / "data" / "contigs"
SUMMARY = Path(__file__).parent / "data" / "summary"

runner = CliRunner()


@pytest.fixture()
def datadir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Fresh tinyamr datadir, published on $GAPIT_DATADIR (MCP resolves it)."""
    target = tmp_path / "datadir"
    shutil.copytree(FIXTURE_DB_DIR, target)
    make_blast_db(target / "tinyamr" / "sequences", "tinyamr")
    monkeypatch.setenv("GAPIT_DATADIR", str(target))
    return target


def exchange(*lines_or_messages: object) -> list[dict[str, Any]]:
    """Feed the loop one line per item (str = raw line, dict = JSON-RPC
    message); return the parsed response objects (parsed-JSON boundary:
    values stay Any until asserted)."""
    raw = "".join(
        (item if isinstance(item, str) else json.dumps(item)) + "\n" for item in lines_or_messages
    )
    out = StringIO()
    serve(StringIO(raw), out)
    return [json.loads(line) for line in out.getvalue().splitlines()]


def tool_call(name: str, arguments: dict[str, object], *, id_value: int = 1) -> dict[str, object]:
    return {
        "jsonrpc": "2.0",
        "id": id_value,
        "method": "tools/call",
        "params": {"name": name, "arguments": arguments},
    }


def test_initialize_echoes_requested_protocol_version() -> None:
    """Given initialize with a requested version, When served, Then the
    response echoes it with gapit serverInfo and tools capabilities."""
    (response,) = exchange(
        {
            "jsonrpc": "2.0",
            "id": 1,
            "method": "initialize",
            "params": {"protocolVersion": "2025-03-26"},
        }
    )
    result = response["result"]
    assert result["protocolVersion"] == "2025-03-26"
    assert result["capabilities"] == {"tools": {}}
    assert result["serverInfo"]["name"] == "gapit"
    assert result["serverInfo"]["version"]


def test_initialize_defaults_when_version_absent_or_empty() -> None:
    """Given initialize without params (or an empty version), When served,
    Then protocolVersion falls back to the 2025-06-18 default."""
    responses = exchange(
        {"jsonrpc": "2.0", "id": "a", "method": "initialize"},
        {"jsonrpc": "2.0", "id": "b", "method": "initialize", "params": {"protocolVersion": ""}},
    )
    assert [r["result"]["protocolVersion"] for r in responses] == [
        "2025-06-18",
        "2025-06-18",
    ]


def test_notifications_get_no_response() -> None:
    """Given requests interleaved with notifications, When served, Then only
    the requests produce response lines."""
    responses = exchange(
        {"jsonrpc": "2.0", "id": 1, "method": "initialize"},
        {"jsonrpc": "2.0", "method": "notifications/initialized"},
        {"jsonrpc": "2.0", "method": "notifications/cancelled"},
    )
    assert len(responses) == 1
    assert responses[0]["id"] == 1


def test_tools_list_advertises_nine_tools_with_schemas() -> None:
    """Given tools/list, When served, Then exactly screen/screen_reads/
    summary/schema/db_list/db_fetch/db_build/db_search/db_outdated, each with
    an object inputSchema carrying properties + required."""
    (response,) = exchange({"jsonrpc": "2.0", "id": 2, "method": "tools/list"})
    tools = response["result"]["tools"]
    assert sorted(tool["name"] for tool in tools) == [
        "db_build",
        "db_fetch",
        "db_list",
        "db_outdated",
        "db_search",
        "schema",
        "screen",
        "screen_reads",
        "summary",
    ]
    by_name = {tool["name"]: tool for tool in tools}
    for tool in by_name.values():
        assert tool["description"]
        assert tool["inputSchema"]["type"] == "object"
        assert "properties" in tool["inputSchema"]
        assert "required" in tool["inputSchema"]
    assert by_name["screen"]["inputSchema"]["required"] == ["files"]
    assert by_name["summary"]["inputSchema"]["required"] == ["files"]
    assert by_name["schema"]["inputSchema"]["required"] == ["name"]
    assert by_name["db_build"]["inputSchema"]["required"] == ["name", "fasta"]
    assert by_name["db_search"]["inputSchema"]["required"] == ["term"]


def test_schema_call_returns_report_schema_text() -> None:
    """Given tools/call schema(report), When served, Then non-error text
    carrying the gapit.report/1 marker and a properties block."""
    (response,) = exchange(tool_call("schema", {"name": "report"}))
    assert response["result"]["isError"] is False
    text = response["result"]["content"][0]["text"]
    assert response["result"]["content"][0]["type"] == "text"
    assert "gapit.report/1" in text
    assert '"properties"' in text


def test_schema_call_unknown_name_is_usage_error() -> None:
    (response,) = exchange(tool_call("schema", {"name": "nope"}))
    assert response["result"]["isError"] is True
    assert json.loads(response["result"]["content"][0]["text"])["code"] == "USAGE_ERROR"


def test_screen_call_returns_report_json_with_expected_hit(datadir: Path) -> None:
    """Given tools/call screen on the tinyamr fixture, When served, Then
    non-error text parsing as gapit.report/1 with the tetA hit."""
    (response,) = exchange(
        tool_call("screen", {"files": [str(CONTIGS / "full.fa")], "db": "tinyamr"})
    )
    assert response["result"]["isError"] is False
    document = json.loads(response["result"]["content"][0]["text"])
    assert document["schema"] == "gapit.report/1"
    (hit,) = document["files"][0]["hits"]
    assert hit["gene"] == "tetA"
    assert document["params"]["db"] == "tinyamr"


def test_screen_call_missing_file_is_error_envelope(datadir: Path) -> None:
    """Given tools/call screen with a missing file, When served, Then
    isError=true and the gapit.error/1 envelope code INPUT_NOT_FOUND."""
    (response,) = exchange(
        tool_call("screen", {"files": [str(datadir / "nope.fa")], "db": "tinyamr"})
    )
    assert response["result"]["isError"] is True
    envelope = json.loads(response["result"]["content"][0]["text"])
    assert envelope["schema"] == "gapit.error/1"
    assert envelope["code"] == "INPUT_NOT_FOUND"


def test_screen_call_rejects_missing_files_argument(datadir: Path) -> None:
    (response,) = exchange(tool_call("screen", {"db": "tinyamr"}))
    assert response["result"]["isError"] is True
    assert json.loads(response["result"]["content"][0]["text"])["code"] == "USAGE_ERROR"


def test_screen_call_honors_datadir_argument_without_env(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given NO $GAPIT_DATADIR in the environment, When tools/call screen
    passes the datadir argument, Then the fixture database screens and the
    tetA hit returns — the db tools' datadir symmetry, no env bridge."""
    monkeypatch.delenv("GAPIT_DATADIR", raising=False)
    target = tmp_path / "datadir"
    shutil.copytree(FIXTURE_DB_DIR, target)
    make_blast_db(target / "tinyamr" / "sequences", "tinyamr")
    (response,) = exchange(
        tool_call(
            "screen", {"files": [str(CONTIGS / "full.fa")], "db": "tinyamr", "datadir": str(target)}
        )
    )
    assert response["result"]["isError"] is False
    document = json.loads(response["result"]["content"][0]["text"])
    assert document["schema"] == "gapit.report/1"
    (hit,) = document["files"][0]["hits"]
    assert hit["gene"] == "tetA"


def test_summary_call_returns_summary_json() -> None:
    (response,) = exchange(
        tool_call(
            "summary", {"files": [str(SUMMARY / "sample_a.tsv"), str(SUMMARY / "sample_b.tsv")]}
        )
    )
    assert response["result"]["isError"] is False
    document = json.loads(response["result"]["content"][0]["text"])
    assert document["schema"] == "gapit.summary/1"
    assert len(document["rows"]) == 2


def test_db_list_call_lists_providers(datadir: Path) -> None:
    (response,) = exchange(tool_call("db_list", {}))
    assert response["result"]["isError"] is False
    document = json.loads(response["result"]["content"][0]["text"])
    assert document["schema"] == "gapit.dblist/1"
    names = [provider["name"] for provider in document["providers"]]
    assert "card" in names
    assert "ncbi" in names


def test_unknown_method_returns_32601() -> None:
    (response,) = exchange({"jsonrpc": "2.0", "id": 7, "method": "resources/list"})
    assert response["id"] == 7
    assert response["error"]["code"] == -32601


def test_unknown_tool_returns_32602() -> None:
    (response,) = exchange(tool_call("nope", {}))
    assert response["error"]["code"] == -32602


def test_malformed_line_is_ignored_and_loop_continues() -> None:
    """Given a non-JSON line before a valid request, When served, Then the
    junk is skipped silently and the request still gets its response."""
    responses = exchange(
        "this is not json",
        {"jsonrpc": "2.0", "id": 3, "method": "tools/list"},
    )
    assert len(responses) == 1
    assert len(responses[0]["result"]["tools"]) == 9


def test_eof_terminates_cleanly() -> None:
    """Given an immediately-closed stdin, When served, Then a clean return
    with no output."""
    out = StringIO()
    serve(StringIO(""), out)
    assert out.getvalue() == ""


def test_cli_mcp_subcommand_streams_protocol_lines() -> None:
    """Given `gapit mcp` fed two request lines via CliRunner, When invoked,
    Then exit 0 and exactly two JSON-RPC response lines on stdout."""
    lines = (
        json.dumps({"jsonrpc": "2.0", "id": 1, "method": "initialize"})
        + "\n"
        + json.dumps({"jsonrpc": "2.0", "id": 2, "method": "tools/list"})
        + "\n"
    )
    result = runner.invoke(app, ["mcp"], input=lines)
    assert result.exit_code == 0
    responses = [json.loads(line) for line in result.stdout.splitlines()]
    assert [response["id"] for response in responses] == [1, 2]
    assert responses[0]["result"]["serverInfo"]["name"] == "gapit"
    assert len(responses[1]["result"]["tools"]) == 9
