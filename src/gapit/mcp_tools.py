"""MCP tool implementations (protocol: gapit.mcp; declarations: gapit.mcp_schemas).

Every tool mirrors its CLI twin by calling the SAME shared callables —
never a reimplementation. Failures raise typed errors; the protocol layer
renders them as isError envelopes.
"""

import json
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from gapit import config
from gapit.db_build_ops import perform_build
from gapit.db_ops import db_list_json, perform_fetch
from gapit.db_query_ops import (
    DEFAULT_LIMIT,
    DEFAULT_STALE_DAYS,
    SearchField,
    outdated_tsv_lines,
    perform_outdated,
    perform_search,
)
from gapit.errors import usage_fail
from gapit.formats.schemas import SCHEMA_MODELS
from gapit.formats.summary import render_summary_json
from gapit.reads import ReadTypeEnum
from gapit.screening import AlignerEnum, OutputFormat, run_screen
from gapit.screening_reads import run_screen_assemblies, run_screen_reads
from gapit.summary import SummaryParams, build_summary

# Parsed JSON-RPC argument containers stay dict[str, Any] until the
# per-field isinstance checks in the helpers below pin concrete types (the
# sanctioned Any boundary, card.py precedent).


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


def _optional_string(arguments: dict[str, Any], key: str) -> str | None:
    value = arguments.get(key)
    if value is not None and not isinstance(value, str):
        usage_fail(f"{key} must be a string")
    return value


def _required_string(arguments: dict[str, Any], key: str) -> str:
    value = arguments.get(key)
    if not isinstance(value, str) or not value:
        usage_fail(f"{key} must be a non-empty string")
    return value


def _optional_path(arguments: dict[str, Any], key: str) -> Path | None:
    value = _optional_string(arguments, key)
    return None if value is None else Path(value)


def _integer(arguments: dict[str, Any], key: str, default: int) -> int:
    value = arguments.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        usage_fail(f"{key} must be a non-negative integer")
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


def _optional_output_format(name: str) -> OutputFormat | None:
    """Map an already-validated format name onto the reads use-case's
    OutputFormat|None (None = the json default, SPEC.md §10)."""
    return None if name == "json" else OutputFormat(name)


def _tool_screen(arguments: dict[str, Any]) -> str:
    files = _paths(arguments, "files")
    if not files:
        usage_fail("no input files given (files is required)")
    db_name = _string(arguments, "db", "ncbi")
    datadir = _optional_path(arguments, "datadir")
    minid = _number(arguments, "minid", 80.0)
    mincov = _number(arguments, "mincov", 80.0)
    output_format = _string(arguments, "format", "json")
    if output_format not in ("json", "tsv", "md"):
        usage_fail(f"format must be one of json, tsv, md: got {output_format}")
    aligner_name = _string(arguments, "aligner", "blastn")
    try:
        aligner = AlignerEnum(aligner_name)
    except ValueError:
        usage_fail(f"aligner must be blastn or minimap2: got {aligner_name}")
    min_breadth = _number(arguments, "min_breadth", 90.0)
    min_identity = _number(arguments, "min_identity", 0.0)
    min_mapq = _integer(arguments, "min_mapq", 0)
    if aligner is AlignerEnum.minimap2:
        # tsv is rejected by the use-case itself (reads-mode format rule);
        # non-default minid/mincov too (blastn-only thresholds).
        return run_screen_assemblies(
            files,
            None,
            db_name,
            datadir,
            None,
            min_breadth,
            min_identity,
            min_mapq,
            threads=1,
            jobs=1,
            noheader=False,
            nopath=False,
            output_format=_optional_output_format(output_format),
            quiet=True,
            minid=minid,
            mincov=mincov,
        )
    if (min_breadth, min_identity, min_mapq) != (90.0, 0.0, 0):
        usage_fail("reads-mode parameters require aligner minimap2")
    return run_screen(
        files=files,
        db_name=db_name,
        datadir=datadir,
        minid=minid,
        mincov=mincov,
        threads=1,
        jobs=1,
        fofn=None,
        quiet=True,
        noheader=False,
        nopath=False,
        debug=False,
        output_format=OutputFormat(output_format),
    )


def _reads_paths(arguments: dict[str, Any], key: str) -> list[Path]:
    """Reads path array -> the native list the use-case consumes (no
    join/split round trip, so commas in filenames survive). Empty elements
    are usage errors: an empty string would silently become the cwd."""
    _paths(arguments, key)
    # list[str] proven by the _paths isinstance checks above
    raw: list[str] = arguments.get(key, [])
    if any(item == "" for item in raw):
        usage_fail(f"{key} contains an empty element")
    return [Path(item) for item in raw]


def _tool_screen_reads(arguments: dict[str, Any]) -> str:
    r1 = _reads_paths(arguments, "r1")
    if not r1:
        usage_fail("no reads files given (r1 is required)")
    r2 = _reads_paths(arguments, "r2") or None
    read_type_name = _string(arguments, "read_type", ReadTypeEnum.sr.value)
    try:
        read_type = ReadTypeEnum(read_type_name)
    except ValueError:
        usage_fail(f"read_type must be one of {', '.join(t.value for t in ReadTypeEnum)}")
    output_format = _string(arguments, "format", "json")
    if output_format not in ("json", "md"):
        usage_fail(f"format must be one of json, md: got {output_format}")
    min_breadth = _number(arguments, "min_breadth", 90.0)
    if not 0.0 <= min_breadth <= 100.0:
        usage_fail(f"min_breadth must be in [0, 100]: got {min_breadth}")
    min_identity = _number(arguments, "min_identity", 0.0)
    if not 0.0 <= min_identity <= 100.0:
        usage_fail(f"min_identity must be in [0, 100]: got {min_identity}")
    min_mapq = _integer(arguments, "min_mapq", 0)
    # lane pairing is validated by the use-case (frozen CLI message)
    return run_screen_reads(
        r1,
        r2,
        _string(arguments, "db", "ncbi"),
        _optional_path(arguments, "datadir"),
        read_type,
        min_breadth,
        min_identity,
        min_mapq,
        threads=1,
        output_format=_optional_output_format(output_format),
        quiet=True,
    )


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
    model = SCHEMA_MODELS.get(name)
    if model is None:
        usage_fail(f"unknown schema name: {name} (choose from: {', '.join(SCHEMA_MODELS)})")
    return json.dumps(model.model_json_schema(by_alias=True), indent=2)


def _tool_db_list(arguments: dict[str, Any]) -> str:
    return db_list_json(config.resolve_datadir(None))


def _tool_db_fetch(arguments: dict[str, Any]) -> str:
    receipts = perform_fetch(
        _optional_string(arguments, "name"),
        _optional_path(arguments, "datadir"),
        force=_flag(arguments, "force"),
    )
    return "\n".join(receipt.model_dump_json() for receipt in receipts)


def _tool_db_build(arguments: dict[str, Any]) -> str:
    dbtype = _optional_string(arguments, "dbtype")
    if dbtype is not None and dbtype not in ("nucl", "prot"):
        usage_fail(f"dbtype must be nucl or prot: got {dbtype}")
    receipt = perform_build(
        _required_string(arguments, "name"),
        Path(_required_string(arguments, "fasta")),
        _optional_path(arguments, "tsv"),
        dbtype,
        _string(arguments, "description", ""),
        _optional_path(arguments, "datadir"),
        _flag(arguments, "force"),
        warn=lambda message: None,
    )
    return receipt.model_dump_json()


def _tool_db_search(arguments: dict[str, Any]) -> str:
    try:
        field = SearchField(_string(arguments, "field", "any"))
    except ValueError:
        usage_fail(f"field must be one of {', '.join(f.value for f in SearchField)}")
    hits, _total = perform_search(
        _required_string(arguments, "term"),
        _optional_path(arguments, "datadir"),
        db=_optional_string(arguments, "db"),
        field=field,
        exact=_flag(arguments, "exact"),
        limit=_integer(arguments, "limit", DEFAULT_LIMIT),
    )
    return "\n".join(hits)


def _tool_db_outdated(arguments: dict[str, Any]) -> str:
    days = _integer(arguments, "days", DEFAULT_STALE_DAYS)
    lines = outdated_tsv_lines(perform_outdated(_optional_path(arguments, "datadir"), days=days))
    return "\n".join(lines)


TOOL_HANDLERS: dict[str, Callable[[dict[str, Any]], str]] = {
    "screen": _tool_screen,
    "screen_reads": _tool_screen_reads,
    "summary": _tool_summary,
    "schema": _tool_schema,
    "db_list": _tool_db_list,
    "db_fetch": _tool_db_fetch,
    "db_build": _tool_db_build,
    "db_search": _tool_db_search,
    "db_outdated": _tool_db_outdated,
}
