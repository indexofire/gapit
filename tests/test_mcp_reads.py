"""MCP reads tools: screen_reads (FASTQ via minimap2) and screen's
aligner=minimap2 assembly survey, driven through the real serve loop.

Same offline fixtures as the CLI reads suites (tests/test_cli_reads.py):
tinyreads needs no BLAST index; the tinyamr survey fixture mirrors
test_mcp.py (makeblastdb). minimap2 soft-clips the short tinyamr genes
(tests/test_reads_integration.py), so the survey asserts the document
contract, not gene presence.
"""

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

READS_DB_DIR = Path(__file__).parent / "data" / "reads_db"
AMR_DB_DIR = Path(__file__).parent / "data" / "db"
READS = Path(__file__).parent / "data" / "reads"
CONTIGS = Path(__file__).parent / "data" / "contigs"

runner = CliRunner()


def exchange(*lines_or_messages: object) -> list[dict[str, Any]]:
    """Feed the loop one line per item (str = raw line, dict = JSON-RPC
    message); return the parsed response objects."""
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


def call_tool(name: str, arguments: dict[str, object]) -> tuple[bool, str]:
    """One tools/call round trip -> (isError, text)."""
    (response,) = exchange(tool_call(name, arguments))
    result = response["result"]
    return result["isError"], result["content"][0]["text"]


def envelope_code(text: str) -> str:
    return json.loads(text)["code"]


@pytest.fixture()
def reads_datadir(tmp_path: Path) -> Path:
    """Fresh tinyreads datadir (minimap2 rows need no BLAST index)."""
    target = tmp_path / "datadir"
    shutil.copytree(READS_DB_DIR, target)
    return target


@pytest.fixture()
def amr_datadir(tmp_path: Path) -> Path:
    """Fresh tinyamr datadir with a BLAST index (test_mcp.py recipe)."""
    target = tmp_path / "datadir"
    shutil.copytree(AMR_DB_DIR, target)
    make_blast_db(target / "tinyamr" / "sequences", "tinyamr")
    return target


def reads_args(reads_datadir: Path, **overrides: object) -> dict[str, object]:
    """Default screen_reads arguments for one paired tetX lane."""
    args: dict[str, object] = {
        "r1": [str(READS / "tetx_R1.fq")],
        "r2": [str(READS / "tetx_R2.fq")],
        "db": "tinyreads",
        "datadir": str(reads_datadir),
    }
    args.update(overrides)
    return args


# ----------------------------------------------------------- screen_reads --


def test_screen_reads_paired_fastq_returns_reads1_json(reads_datadir: Path) -> None:
    """Given tools/call screen_reads on one paired tetX lane, When served,
    Then non-error text parsing as gapit.reads/1 with the sr preset and the
    tetX gene present (unmapped db genes are omitted from the document)."""
    is_error, text = call_tool("screen_reads", reads_args(reads_datadir))
    assert is_error is False
    document = json.loads(text)
    assert document["schema"] == "gapit.reads/1"
    assert document["params"]["db"] == "tinyreads"
    assert document["params"]["read_type"] == "sr"
    genes = {entry["gene"]: entry["present"] for entry in document["files"][0]["genes"]}
    assert genes == {"tetX": True}


def test_screen_reads_min_identity_emits_reads2(reads_datadir: Path) -> None:
    """Given screen_reads with min_identity > 0, When served, Then the
    document is the filtered gapit.reads/2 schema."""
    is_error, text = call_tool("screen_reads", reads_args(reads_datadir, min_identity=90.0))
    assert is_error is False
    document = json.loads(text)
    assert document["schema"] == "gapit.reads/2"
    assert document["params"]["min_identity"] == 90.0


def test_screen_reads_md_format_renders_markdown(reads_datadir: Path) -> None:
    """Given screen_reads with format md, When served, Then non-error text
    carrying the gapit.reads/1 YAML frontmatter."""
    is_error, text = call_tool("screen_reads", reads_args(reads_datadir, format="md"))
    assert is_error is False
    assert text.startswith("---\n")
    assert "schema: gapit.reads/1" in text


def test_screen_reads_rejects_tsv_format(reads_datadir: Path) -> None:
    """Given screen_reads with format tsv (CLI-only), When served, Then
    isError with the USAGE_ERROR envelope."""
    is_error, text = call_tool("screen_reads", reads_args(reads_datadir, format="tsv"))
    assert is_error is True
    assert envelope_code(text) == "USAGE_ERROR"


def test_screen_reads_r2_lane_mismatch_is_usage_error(reads_datadir: Path) -> None:
    """Given two r1 lanes but one r2 mate, When served, Then isError with
    the frozen lanes-must-pair-up message."""
    is_error, text = call_tool(
        "screen_reads",
        reads_args(
            reads_datadir,
            r1=[str(READS / "tetx_lane1.fq"), str(READS / "tetx_lane2.fq")],
            r2=[str(READS / "tetx_R2.fq")],
        ),
    )
    assert is_error is True
    assert envelope_code(text) == "USAGE_ERROR"
    assert "lanes must pair up" in text


def test_screen_reads_invalid_read_type_is_usage_error(reads_datadir: Path) -> None:
    is_error, text = call_tool("screen_reads", reads_args(reads_datadir, read_type="pacbio"))
    assert is_error is True
    assert envelope_code(text) == "USAGE_ERROR"


def test_screen_reads_empty_r1_is_usage_error(reads_datadir: Path) -> None:
    is_error, text = call_tool("screen_reads", reads_args(reads_datadir, r1=[]))
    assert is_error is True
    assert envelope_code(text) == "USAGE_ERROR"


def test_screen_reads_empty_path_element_is_usage_error(reads_datadir: Path) -> None:
    """Given r1 containing an empty path string, When served, Then isError
    (an empty element would silently become the working directory)."""
    is_error, text = call_tool("screen_reads", reads_args(reads_datadir, r1=[""]))
    assert is_error is True
    assert envelope_code(text) == "USAGE_ERROR"


# -------------------------------------------------- screen aligner matrix --


def test_screen_aligner_minimap2_surveys_assembly(amr_datadir: Path) -> None:
    """Given tools/call screen with aligner minimap2 on a contigs file, When
    served, Then the reads-engine survey returns gapit.reads/1 with the
    auto-resolved map-ont preset."""
    is_error, text = call_tool(
        "screen",
        {
            "files": [str(CONTIGS / "full.fa")],
            "db": "tinyamr",
            "datadir": str(amr_datadir),
            "aligner": "minimap2",
        },
    )
    assert is_error is False
    document = json.loads(text)
    assert document["schema"] == "gapit.reads/1"
    assert document["params"]["read_type"] == "map-ont"


def test_screen_aligner_blastn_rejects_reads_parameters(amr_datadir: Path) -> None:
    """Given screen without aligner (blastn default) but with a reads-mode
    threshold, When served, Then isError naming the aligner requirement."""
    is_error, text = call_tool(
        "screen",
        {
            "files": [str(CONTIGS / "full.fa")],
            "db": "tinyamr",
            "datadir": str(amr_datadir),
            "min_identity": 90.0,
        },
    )
    assert is_error is True
    assert envelope_code(text) == "USAGE_ERROR"
    assert "reads-mode parameters require aligner minimap2" in text


# -------------------------------------------------------------- cli smoke --


def test_cli_mcp_serve_streams_screen_reads_response(reads_datadir: Path) -> None:
    """Given `gapit mcp` fed initialize + a single-end screen_reads call via
    CliRunner, When invoked, Then exit 0 and the gapit.reads/1 document
    rides the second response line (r2 is optional)."""
    request = tool_call(
        "screen_reads",
        {
            "r1": [str(READS / "tetx_full.fq")],
            "db": "tinyreads",
            "datadir": str(reads_datadir),
        },
        id_value=2,
    )
    lines = (
        json.dumps({"jsonrpc": "2.0", "id": 1, "method": "initialize"})
        + "\n"
        + json.dumps(request)
        + "\n"
    )
    result = runner.invoke(app, ["mcp"], input=lines)
    assert result.exit_code == 0
    responses = [json.loads(line) for line in result.stdout.splitlines()]
    assert [response["id"] for response in responses] == [1, 2]
    assert responses[1]["result"]["isError"] is False
    document = json.loads(responses[1]["result"]["content"][0]["text"])
    assert document["schema"] == "gapit.reads/1"
    genes = {entry["gene"]: entry["present"] for entry in document["files"][0]["genes"]}
    assert genes["tetX"] is True
