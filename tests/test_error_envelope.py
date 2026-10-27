"""Error envelope tests: every typed failure renders gapit.error/1 on stderr."""

import ast
import json
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner
from pydantic import TypeAdapter

from gapit.cli import app
from gapit.db import make_blast_db
from gapit.errors import ERROR_CODES, ErrorEnvelope, InputError, ensure_input_file, render_error

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"
CONTIGS = Path(__file__).parent / "data" / "contigs"

runner = CliRunner()
envelope_adapter = TypeAdapter(ErrorEnvelope)


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    """A datadir holding an indexed copy of tinyamr."""
    target = tmp_path / "datadir" / "tinyamr"
    target.mkdir(parents=True)
    shutil.copy(FIXTURE_DB_DIR / "tinyamr" / "sequences", target / "sequences")
    make_blast_db(target / "sequences", "tinyamr")
    return target.parent


def last_envelope(stderr: str) -> ErrorEnvelope:
    """Parse the last non-empty stderr line as the error envelope."""
    lines = [line for line in stderr.splitlines() if line.strip()]
    assert lines, "expected an error envelope on stderr"
    return envelope_adapter.validate_json(lines[-1])


def test_envelope_is_single_line_with_four_keys(tmp_path: Path) -> None:
    """Given any typed failure, When rendered, Then stderr's last line is a
    single-line JSON object with exactly the four documented keys in order."""
    result = runner.invoke(app, ["db", "outdated", "--datadir", str(tmp_path / "nope")])
    assert result.exit_code == 4
    lines = [line for line in result.stderr.splitlines() if line.strip()]
    assert len(lines) == 1
    assert lines[0].startswith('{"schema":"gapit.error/1","code":')
    payload = json.loads(lines[0])
    assert list(payload) == ["schema", "code", "message", "context"]


def test_missing_datadir_envelope(tmp_path: Path) -> None:
    result = runner.invoke(app, ["db", "outdated", "--datadir", str(tmp_path / "nope")])
    assert result.exit_code == 4
    envelope = last_envelope(result.stderr)
    assert envelope.code == "DATADIR_NOT_FOUND"
    assert envelope.context["datadir"] == str(tmp_path / "nope")


def test_unknown_db_envelope(datadir: Path) -> None:
    result = runner.invoke(
        app,
        ["screen", "--db", "nope", "--datadir", str(datadir), str(CONTIGS / "full.fa")],
    )
    assert result.exit_code == 4
    envelope = last_envelope(result.stderr)
    assert envelope.code == "DATABASE_NOT_FOUND"
    assert envelope.context["db"] == "nope"
    assert "tinyamr" in envelope.message


def test_missing_input_file_envelope(datadir: Path) -> None:
    result = runner.invoke(
        app,
        ["screen", "--db", "tinyamr", "--datadir", str(datadir), str(datadir / "nope.fa")],
    )
    assert result.exit_code == 5
    envelope = last_envelope(result.stderr)
    assert envelope.code == "INPUT_NOT_FOUND"
    assert envelope.context["file"] == str(datadir / "nope.fa")


def test_ensure_input_file_parametrizes_the_kind(tmp_path: Path) -> None:
    """Given the shared existence helper, When a missing path is checked with
    a kind, Then InputError carries that kind in the message plus the
    INPUT_NOT_FOUND code and file context (bytes match the former inline
    raises in screening.py / screening_reads.py)."""
    missing = tmp_path / "nope.fq"
    with pytest.raises(InputError) as excinfo:
        ensure_input_file(missing, "reads file")
    assert str(excinfo.value) == f"reads file not found or unreadable: {missing}"
    assert excinfo.value.code == "INPUT_NOT_FOUND"
    assert excinfo.value.context == {"file": str(missing)}
    with pytest.raises(InputError) as default:
        ensure_input_file(missing)
    assert str(default.value) == f"input file not found or unreadable: {missing}"


def test_junk_input_envelope(datadir: Path, tmp_path: Path) -> None:
    junk = tmp_path / "junk.txt"
    junk.write_text("this is not sequence data at all\n", encoding="utf-8")
    result = runner.invoke(app, ["screen", "--db", "tinyamr", "--datadir", str(datadir), str(junk)])
    assert result.exit_code == 5
    envelope = last_envelope(result.stderr)
    assert envelope.code == "INVALID_INPUT"


def test_invalid_minid_usage_envelope(datadir: Path) -> None:
    result = runner.invoke(
        app,
        [
            "screen",
            "--db",
            "tinyamr",
            "--minid",
            "0",
            "--datadir",
            str(datadir),
            str(CONTIGS / "full.fa"),
        ],
    )
    assert result.exit_code == 2
    envelope = last_envelope(result.stderr)
    assert envelope.code == "USAGE_ERROR"
    assert "--minid" in envelope.message


def test_unexpected_exception_envelope() -> None:
    """Given a non-GapitError failure, When rendered, Then the UNEXPECTED
    fallback envelope with exit 1 is used."""
    line = render_error(ValueError("boom"))
    envelope = envelope_adapter.validate_json(line)
    assert envelope.code == "UNEXPECTED"
    assert "ValueError" in envelope.message
    assert "boom" in envelope.message
    assert envelope.context == {}


def test_unknown_option_renders_usage_envelope() -> None:
    """Given an unknown option, When click rejects it before the command body,
    Then stderr is still the gapit.error/1 usage envelope, not click prose."""
    result = runner.invoke(app, ["screen", "--bogus"])
    assert result.exit_code == 2
    envelope = last_envelope(result.stderr)
    assert envelope.code == "USAGE_ERROR"
    assert "--bogus" in envelope.message


def test_missing_required_option_renders_usage_envelope() -> None:
    """Given a missing required --db, When click rejects the invocation,
    Then the same typed usage envelope is emitted (exit 2)."""
    result = runner.invoke(app, ["screen", str(CONTIGS / "full.fa")])
    assert result.exit_code == 2
    envelope = last_envelope(result.stderr)
    assert envelope.code == "USAGE_ERROR"
    assert "--db" in envelope.message


def test_error_schema_exposes_the_registered_code_enum() -> None:
    """Given gapit.error/1, When introspected, Then code carries the full
    machine-readable error-code enum (no doc-only list)."""
    schema = ErrorEnvelope.model_json_schema(by_alias=True)
    assert set(schema["properties"]["code"]["enum"]) == set(ERROR_CODES)


def test_every_raise_site_uses_a_registered_code() -> None:
    """Given every typed raise helper/call in src/gapit, When scanned, Then
    its code is in ERROR_CODES (so the public enum cannot drift)."""
    source = Path(__file__).parents[1] / "src" / "gapit"
    names = {"_fail", "GapitError", "DatabaseError", "InputError", "UsageError", "DependencyError"}
    found: set[str] = set()
    for path in source.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            function = node.func
            name = (
                function.id
                if isinstance(function, ast.Name)
                else function.attr
                if isinstance(function, ast.Attribute)
                else ""
            )
            if name in names and node.args and isinstance(node.args[0], ast.Constant):
                value = node.args[0].value
                if isinstance(value, str) and value.isupper():
                    found.add(value)
            for keyword in node.keywords:
                if keyword.arg != "code" or not isinstance(keyword.value, ast.Constant):
                    continue
                value = keyword.value.value
                if isinstance(value, str):
                    found.add(value)
    assert found <= ERROR_CODES, f"unregistered error codes: {sorted(found - ERROR_CODES)}"
