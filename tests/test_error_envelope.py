"""Error envelope tests: every typed failure renders gapit.error/1 on stderr."""

import json
import shutil
from pathlib import Path

import pytest
from pydantic import TypeAdapter
from typer.testing import CliRunner

from gapit.cli import app
from gapit.db import make_blast_db
from gapit.errors import ErrorEnvelope, render_error

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
    result = runner.invoke(app, ["list", "--datadir", str(tmp_path / "nope")])
    assert result.exit_code == 4
    lines = [line for line in result.stderr.splitlines() if line.strip()]
    assert len(lines) == 1
    assert lines[0].startswith('{"schema":"gapit.error/1","code":')
    payload = json.loads(lines[0])
    assert list(payload) == ["schema", "code", "message", "context"]


def test_missing_datadir_envelope(tmp_path: Path) -> None:
    result = runner.invoke(app, ["list", "--datadir", str(tmp_path / "nope")])
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
    result = runner.invoke(app, ["screen", "--datadir", str(datadir), str(datadir / "nope.fa")])
    assert result.exit_code == 5
    envelope = last_envelope(result.stderr)
    assert envelope.code == "INPUT_NOT_FOUND"
    assert envelope.context["file"] == str(datadir / "nope.fa")


def test_junk_input_envelope(datadir: Path, tmp_path: Path) -> None:
    junk = tmp_path / "junk.txt"
    junk.write_text("this is not sequence data at all\n", encoding="utf-8")
    result = runner.invoke(app, ["screen", "--db", "tinyamr", "--datadir", str(datadir), str(junk)])
    assert result.exit_code == 5
    envelope = last_envelope(result.stderr)
    assert envelope.code == "INVALID_INPUT"


def test_invalid_minid_usage_envelope(datadir: Path) -> None:
    result = runner.invoke(
        app, ["screen", "--minid", "0", "--datadir", str(datadir), str(CONTIGS / "full.fa")]
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
