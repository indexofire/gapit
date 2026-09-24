r"""Minimal MCP stdio server: gapit as a tool provider for agents.

Hand-rolled on purpose (AGENTS.md §2 keeps the dependency list deliberately
short): no `mcp` SDK — just the essential MCP stdio behavior, JSON-RPC 2.0,
one message per line on stdin and one response line on stdout. Limitations:
single messages only (no batch arrays); non-JSON lines are ignored silently
(robustness over -32700); JSON lines that fail frame validation get a
-32600/-32602 error when they carry a request id, so no client ever hangs.
This module is the PROTOCOL only — frames,
dispatch, and the serve loop; the tool implementations live in
:mod:`gapit.mcp_tools` and the tools/list declarations (inputSchemas) in
:mod:`gapit.mcp_schemas` (screen, summary, schema, db_list, db_fetch,
db_build, db_search, db_outdated). Tool failures return isError=true with
the gapit.error/1 envelope as text; stdout is protocol-only.
"""

import json
import sys
from collections.abc import Iterable
from typing import Any, TextIO, TypeGuard

import typer
from pydantic import BaseModel, ConfigDict, ValidationError

from gapit import __version__
from gapit.errors import render_error
from gapit.mcp_schemas import TOOLS
from gapit.mcp_tools import TOOL_HANDLERS

# Parsed JSON-RPC values are the one sanctioned Any boundary (card.py
# precedent): frame/argument containers stay dict[str, Any] until the
# per-field checks in mcp_tools pin concrete types.
JsonRpcId = str | int | float | None

PROTOCOL_VERSION = "2025-06-18"


class _Frame(BaseModel, frozen=True):
    """One stdin frame (request or notification); ``id`` presence in
    model_fields_set distinguishes requests, per JSON-RPC 2.0."""

    model_config = ConfigDict(extra="ignore")

    method: str | None = None
    params: dict[str, Any] | None = None
    id: JsonRpcId = None


class _ToolCall(BaseModel, frozen=True):
    """tools/call params: the tool ``name`` and its ``arguments`` object."""

    model_config = ConfigDict(extra="ignore")

    name: str = ""
    arguments: dict[str, Any] | None = None


def _result(id_value: JsonRpcId, result: dict[str, object]) -> dict[str, object]:
    return {"jsonrpc": "2.0", "id": id_value, "result": result}


def _error(id_value: JsonRpcId, code: int, message: str) -> dict[str, object]:
    return {"jsonrpc": "2.0", "id": id_value, "error": {"code": code, "message": message}}


def _text(text: str, *, is_error: bool) -> dict[str, object]:
    return {"content": [{"type": "text", "text": text}], "isError": is_error}


def _initialize(params: dict[str, Any] | None) -> dict[str, object]:
    version = (params or {}).get("protocolVersion", "")
    requested = version if isinstance(version, str) else ""
    return {
        "protocolVersion": requested or PROTOCOL_VERSION,
        "capabilities": {"tools": {}},
        "serverInfo": {"name": "gapit", "version": __version__},
    }


def _tools_call(id_value: JsonRpcId, params: dict[str, Any] | None) -> dict[str, object]:
    try:
        call = _ToolCall.model_validate(params or {})
    except ValidationError:
        return _error(id_value, -32602, "invalid tools/call params")
    tool = TOOL_HANDLERS.get(call.name)
    if tool is None:
        return _error(id_value, -32602, f"unknown tool: {call.name or '(none)'}")
    try:
        text = tool(call.arguments or {})
    except Exception as exc:  # gapit.dispatch semantics: envelope, not a crash
        return _result(id_value, _text(render_error(exc), is_error=True))
    return _result(id_value, _text(text, is_error=False))


def _is_json_object(value: Any) -> TypeGuard[dict[str, Any]]:
    """Narrow a parsed JSON value to the sanctioned dict[str, Any] boundary
    (plain isinstance would surface dict[Unknown, Unknown])."""
    return isinstance(value, dict)


def _reject(value: Any) -> dict[str, object] | None:
    """A line that parsed as JSON but failed _Frame validation: answer the
    request id with -32600 (missing/non-string method) or -32602 (params not
    an object); id-less messages stay ignored like notifications."""
    if not _is_json_object(value):
        return None
    id_value = value.get("id")
    if id_value is None:
        return None
    if not isinstance(value.get("method"), str):
        return _error(id_value, -32600, "invalid request: method must be a string")
    return _error(id_value, -32602, "invalid request: params must be an object")


def _handle(frame: _Frame) -> dict[str, object] | None:
    """One parsed frame -> one response object; None = no response."""
    if "id" not in frame.model_fields_set:
        return None  # notifications never get responses
    if frame.method is None:
        return _error(frame.id, -32600, "invalid request: method is required")
    id_value = frame.id
    if frame.method == "initialize":
        return _result(id_value, _initialize(frame.params))
    if frame.method == "tools/list":
        return _result(id_value, {"tools": list(TOOLS)})
    if frame.method == "tools/call":
        return _tools_call(id_value, frame.params)
    return _error(id_value, -32601, f"method not found: {frame.method}")


def serve(stdin: Iterable[str], stdout: TextIO) -> None:
    """Serve newline-delimited JSON-RPC 2.0 until EOF; non-JSON lines and
    id-less invalid messages are ignored, requests that fail frame
    validation get a -32600/-32602 error instead of silence."""
    for line in stdin:
        try:
            value: object = json.loads(line)
        except json.JSONDecodeError:
            continue
        try:
            frame = _Frame.model_validate(value)
        except ValidationError:
            response = _reject(value)
        else:
            response = _handle(frame)
        if response is not None:
            stdout.write(json.dumps(response, separators=(",", ":")) + "\n")
            stdout.flush()


def register_mcp_command(app: typer.Typer) -> None:
    """Attach the `mcp` subcommand (cli.py stays import + registration only)."""
    app.command("mcp")(main)


def main() -> None:
    """`gapit-mcp` console script; also backs the `gapit mcp` subcommand."""
    serve(sys.stdin, sys.stdout)
