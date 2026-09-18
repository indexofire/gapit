r"""Minimal MCP stdio server: gapit as a read-only tool provider for agents.

Hand-rolled on purpose (AGENTS.md §2 keeps the dependency list deliberately
short): no `mcp` SDK — just the essential MCP stdio behavior, JSON-RPC 2.0,
one message per line on stdin and one response line on stdout. Limitations:
single messages only (no batch arrays); non-JSON lines are ignored silently
(robustness over -32700). Tools: screen (gapit.report/1 default), summary
(gapit.summary/1), schema, db_list (gapit.dblist/1). Tool failures return
isError=true with the gapit.error/1 envelope as text; stdout is protocol-only.
"""

import json
import sys
from collections.abc import Callable, Iterable
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, TextIO

import typer
from pydantic import BaseModel, ConfigDict, ValidationError

from gapit import __version__, config, db
from gapit.blast import screen_file
from gapit.cmd_db import DbListDocument, DbListEntry
from gapit.errors import DatabaseError, ErrorEnvelope, InputError, render_error
from gapit.formats import json as fmt_json
from gapit.formats.md import render_markdown
from gapit.formats.summary import SummaryDocument, render_summary_json
from gapit.formats.tsv import format_tsv
from gapit.providers import REGISTRY
from gapit.records import read_manifest
from gapit.report import ScreeningParams
from gapit.screening import usage_fail
from gapit.summary import SummaryParams, build_summary

# Parsed JSON-RPC values are the one sanctioned Any boundary (card.py
# precedent): frame/argument containers stay dict[str, Any] until the
# per-field isinstance checks pin concrete types.
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


# Mirrors cli._SCHEMA_MODELS (private there; basedpyright reportPrivateUsage
# blocks the import — Wave A3 _run local-mirror precedent).
_SCHEMA_MODELS: dict[str, type[BaseModel]] = {
    "report": fmt_json.ReportDocument,
    "reads": fmt_json.ReadsDocument,
    "summary": SummaryDocument,
    "list": fmt_json.ListDocument,
    "error": ErrorEnvelope,
    "version": fmt_json.VersionDocument,
}


def _object_schema(properties: dict[str, object], required: list[str]) -> dict[str, object]:
    return {"type": "object", "properties": properties, "required": required}


_FILES: dict[str, object] = {"type": "array", "items": {"type": "string"}}

_TOOLS: list[dict[str, object]] = [
    {
        "name": "screen",
        "description": "Screen contig files for AMR/virulence genes (json = gapit.report/1).",
        "inputSchema": _object_schema(
            {
                "files": _FILES,
                "db": {"type": "string", "default": "ncbi"},
                "minid": {"type": "number"},
                "mincov": {"type": "number"},
                "format": {"type": "string", "enum": ["json", "tsv", "md"], "default": "json"},
            },
            ["files"],
        ),
    },
    {
        "name": "summary",
        "description": "Summarize report tables into a gapit.summary/1 matrix.",
        "inputSchema": _object_schema(
            {"files": _FILES, "identity": {"type": "boolean"}, "nopath": {"type": "boolean"}},
            ["files"],
        ),
    },
    {
        "name": "schema",
        "description": "Print the JSON Schema of a gapit output document.",
        "inputSchema": _object_schema(
            {"name": {"type": "string", "enum": sorted(_SCHEMA_MODELS)}},
            ["name"],
        ),
    },
    {
        "name": "db_list",
        "description": "List database providers and their installed state (gapit.dblist/1).",
        "inputSchema": _object_schema({}, []),
    },
]


def _paths(arguments: dict[str, Any], key: str) -> list[Path]:
    value: Any = arguments.get(key, [])
    # items stays unnarrowed Any: iterating the isinstance-narrowed value
    # would leak Unknown into basedpyright strict; value proves list-ness.
    items: Any = arguments.get(key, [])
    if not isinstance(value, list) or not all(isinstance(item, str) for item in items):
        usage_fail(f"{key} must be an array of file path strings")
    return [Path(item) for item in items]


def _string(arguments: dict[str, Any], key: str, default: str) -> str:
    value = arguments.get(key, default)
    if not isinstance(value, str):
        usage_fail(f"{key} must be a string")
    return value


def _number(arguments: dict[str, Any], key: str, default: float) -> float:
    value = arguments.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        usage_fail(f"{key} must be a number")
    return float(value)


def _flag(arguments: dict[str, Any], key: str) -> bool:
    value = arguments.get(key, False)
    if not isinstance(value, bool):
        usage_fail(f"{key} must be a boolean")
    return value


def _find_database(datadir: Path, name: str) -> db.Database:
    """Minimal re-derivation of screening._find_database (private there;
    reportPrivateUsage blocks the import) — same errors, same wording."""
    databases = db.discover_databases(datadir)
    for database in databases:
        if database.name == name:
            return database
    available = ", ".join(entry.name for entry in databases) or "(none)"
    raise DatabaseError(
        f"Database {name} is not in {datadir}. Available: {available}",
        code="DATABASE_NOT_FOUND",
        context={"db": name, "datadir": str(datadir)},
    )


def _tool_screen(arguments: dict[str, Any]) -> str:
    files = _paths(arguments, "files")
    if not files:
        usage_fail("no input files given (files is required)")
    db_name = _string(arguments, "db", "ncbi")
    minid = _number(arguments, "minid", 80.0)
    mincov = _number(arguments, "mincov", 80.0)
    if not 0.0 < minid <= 100.0:
        usage_fail(f"minid must be in (0, 100]: got {minid}")
    if not 0.0 <= mincov <= 100.0:
        usage_fail(f"mincov must be in [0, 100]: got {mincov}")
    output_format = _string(arguments, "format", "json")
    if output_format not in ("json", "tsv", "md"):
        usage_fail(f"format must be one of json, tsv, md: got {output_format}")
    for path in files:
        if not path.is_file():
            raise InputError(
                f"input file not found or unreadable: {path}",
                code="INPUT_NOT_FOUND",
                context={"file": str(path)},
            )
    params = ScreeningParams(db=db_name, minid=minid, mincov=mincov, threads=1)
    database = _find_database(config.resolve_datadir(None), db_name)
    reports = [screen_file(path, database, params) for path in files]
    now = datetime.now(UTC)
    if output_format == "json":
        return fmt_json.render_json(reports, params, now=now)
    if output_format == "md":
        return render_markdown(reports, params, now=now)
    return format_tsv(reports, csv=False, noheader=False, nopath=False)


def _tool_summary(arguments: dict[str, Any]) -> str:
    files = _paths(arguments, "files")
    if not files:
        usage_fail("summary needs >= 1 report file(s)")
    params = SummaryParams(identity=_flag(arguments, "identity"), nopath=_flag(arguments, "nopath"))
    # The CLI warns about duplicate inputs on stderr; MCP reserves stderr for
    # protocol-internal errors, so those warnings are dropped here.
    matrix = build_summary(files, params, warn=lambda message: None)
    return render_summary_json(matrix, now=datetime.now(UTC))


def _tool_schema(arguments: dict[str, Any]) -> str:
    name = _string(arguments, "name", "")
    model = _SCHEMA_MODELS.get(name)
    if model is None:
        usage_fail(f"unknown schema name: {name} (choose from: {', '.join(_SCHEMA_MODELS)})")
    return json.dumps(model.model_json_schema(by_alias=True), indent=2)


def _tool_db_list(arguments: dict[str, Any]) -> str:
    """gapit.dblist/1, mirroring cmd_db.db_list_command's JSON branch."""
    root = config.resolve_datadir(None)
    entries: list[DbListEntry] = []
    for name in sorted(REGISTRY):
        manifest = root / name / "gapit-manifest.json"
        installed = manifest.is_file()
        entries.append(
            DbListEntry(
                name=name,
                description=REGISTRY[name].description,
                dbtype=REGISTRY[name].dbtype,
                installed=installed,
                records=read_manifest(manifest).n_records if installed else None,
            )
        )
    return DbListDocument(providers=tuple(entries)).model_dump_json(
        indent=2, by_alias=True, exclude_none=True
    )


_TOOL_HANDLERS: dict[str, Callable[[dict[str, Any]], str]] = {
    "screen": _tool_screen,
    "summary": _tool_summary,
    "schema": _tool_schema,
    "db_list": _tool_db_list,
}


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
    tool = _TOOL_HANDLERS.get(call.name)
    if tool is None:
        return _error(id_value, -32602, f"unknown tool: {call.name or '(none)'}")
    try:
        text = tool(call.arguments or {})
    except Exception as exc:  # CLI _dispatch semantics: envelope, not a crash
        return _result(id_value, _text(render_error(exc), is_error=True))
    return _result(id_value, _text(text, is_error=False))


def _handle(frame: _Frame) -> dict[str, object] | None:
    """One parsed frame -> one response object; None = no response."""
    if frame.method is None or "id" not in frame.model_fields_set:
        return None  # notifications and malformed frames never get responses
    id_value = frame.id
    if frame.method == "initialize":
        return _result(id_value, _initialize(frame.params))
    if frame.method == "tools/list":
        return _result(id_value, {"tools": list(_TOOLS)})
    if frame.method == "tools/call":
        return _tools_call(id_value, frame.params)
    return _error(id_value, -32601, f"method not found: {frame.method}")


def serve(stdin: Iterable[str], stdout: TextIO) -> None:
    """Serve newline-delimited JSON-RPC 2.0 until EOF; non-JSON lines ignored."""
    for line in stdin:
        try:
            frame = _Frame.model_validate_json(line)
        except ValidationError:
            continue
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
