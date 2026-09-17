"""Integration tests: `gapit db fetch` — verified local-file installation.

Generic file plumbing only: plain-text (non-biological) fixtures copied from a
local source path to a target with a streaming SHA256 check. There is no
network fetch and no database/provider semantics here.
"""

import hashlib
import json
from pathlib import Path

import pytest
from pydantic import TypeAdapter
from typer.testing import CliRunner

from gapit.cli import app
from gapit.errors import ErrorEnvelope

runner = CliRunner()
envelope_adapter = TypeAdapter(ErrorEnvelope)

PAYLOAD = b"gapit local install fixture: plain text, not sequence data\n" * 37
PAYLOAD_SHA256 = hashlib.sha256(PAYLOAD).hexdigest()
OLD_TARGET = b"previous install, must survive a failed fetch\n"


def last_envelope(stderr: str) -> ErrorEnvelope:
    """Parse the last non-empty stderr line as the error envelope."""
    lines = [line for line in stderr.splitlines() if line.strip()]
    assert lines, "expected an error envelope on stderr"
    return envelope_adapter.validate_json(lines[-1])


def write_source(tmp_path: Path, name: str = "source.txt") -> Path:
    """A local source file holding the generic fixture payload."""
    source = tmp_path / name
    source.write_bytes(PAYLOAD)
    return source


def test_db_fetch_installs_bytes_for_known_sha256(tmp_path: Path) -> None:
    """Given a local source and its true SHA256, When fetched, Then exit 0,
    the target holds the exact source bytes, stdout is a one-line JSON receipt
    naming the destination and the verified digest, and no temp file remains."""
    source = write_source(tmp_path)
    target_dir = tmp_path / "dest"
    target_dir.mkdir()
    target = target_dir / "installed.txt"
    result = runner.invoke(
        app,
        ["db", "fetch", str(source), "--sha256", PAYLOAD_SHA256, "--output", str(target)],
    )
    assert result.exit_code == 0
    assert target.read_bytes() == PAYLOAD
    receipt = json.loads(result.stdout)
    assert receipt == {"destination": str(target), "sha256": PAYLOAD_SHA256}
    assert [p.name for p in target_dir.iterdir()] == [target.name]


def test_db_fetch_checksum_mismatch_preserves_existing_target(tmp_path: Path) -> None:
    """Given a source whose digest differs from --sha256 and an existing
    target, When fetched, Then exit 5 with a CHECKSUM_MISMATCH envelope on
    stderr, nothing on stdout, the old target bytes intact, and no temp
    leftovers beside it."""
    source = write_source(tmp_path)
    target_dir = tmp_path / "dest"
    target_dir.mkdir()
    target = target_dir / "installed.txt"
    target.write_bytes(OLD_TARGET)
    result = runner.invoke(
        app,
        ["db", "fetch", str(source), "--sha256", "0" * 64, "--output", str(target)],
    )
    assert result.exit_code == 5
    assert result.stdout == ""
    envelope = last_envelope(result.stderr)
    assert envelope.code == "CHECKSUM_MISMATCH"
    assert envelope.context["expected"] == "0" * 64
    assert envelope.context["actual"] == PAYLOAD_SHA256
    assert target.read_bytes() == OLD_TARGET
    assert [p.name for p in target_dir.iterdir()] == [target.name]


@pytest.mark.parametrize("bad", ["deadbeef", "z" * 64, "a" * 63])
def test_db_fetch_rejects_malformed_digest_as_usage_error(tmp_path: Path, bad: str) -> None:
    """Given a --sha256 that is not exactly 64 hex chars, When fetched, Then
    exit 2 with a USAGE_ERROR envelope mentioning --sha256 and no target is
    created."""
    source = write_source(tmp_path)
    target = tmp_path / "installed.txt"
    result = runner.invoke(
        app, ["db", "fetch", str(source), "--sha256", bad, "--output", str(target)]
    )
    assert result.exit_code == 2
    assert result.stdout == ""
    envelope = last_envelope(result.stderr)
    assert envelope.code == "USAGE_ERROR"
    assert "--sha256" in envelope.message
    assert not target.exists()


def test_db_fetch_missing_source_is_input_error(tmp_path: Path) -> None:
    """Given a --output-ready tmpdir but a nonexistent source, When fetched,
    Then exit 5 with an INPUT_NOT_FOUND envelope carrying the source path."""
    source = tmp_path / "absent.txt"
    target = tmp_path / "installed.txt"
    result = runner.invoke(
        app,
        ["db", "fetch", str(source), "--sha256", PAYLOAD_SHA256, "--output", str(target)],
    )
    assert result.exit_code == 5
    assert result.stdout == ""
    envelope = last_envelope(result.stderr)
    assert envelope.code == "INPUT_NOT_FOUND"
    assert envelope.context["source"] == str(source)


def test_db_fetch_accepts_uppercase_digest(tmp_path: Path) -> None:
    """Given the true SHA256 in uppercase, When fetched, Then exit 0 and the
    receipt reports the canonical lowercase digest."""
    source = write_source(tmp_path)
    target = tmp_path / "installed.txt"
    result = runner.invoke(
        app,
        ["db", "fetch", str(source), "--sha256", PAYLOAD_SHA256.upper(), "--output", str(target)],
    )
    assert result.exit_code == 0
    assert json.loads(result.stdout)["sha256"] == PAYLOAD_SHA256
