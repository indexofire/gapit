"""The verified LOCAL-FILE install path: `gapit db install` (from cmd_db.py).

``db install`` checksum-verifies and atomically installs BYTES: SOURCE must
be a path to an existing regular file; no network, no provider IDs, no
archives, no manifests — streaming SHA256 + atomic copy, kept deliberately
independent of db.py and the screening pipeline.

Split from cmd_db.py in Wave G: cmd_db hit the 250 LOC ceiling when fetch
gained the bundled-snapshot defaults; this path is self-contained and no
test imports it directly (everything drives `gapit.cli.app`).
"""

import hashlib
import os
import re
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Annotated

import typer
from pydantic import BaseModel

from gapit.dispatch import dispatch
from gapit.errors import InputError, UsageError

_SHA256_SHAPE = re.compile(r"[0-9a-fA-F]{64}")
_CHUNK_BYTES = 1 << 20


class FetchReceipt(BaseModel, frozen=True):
    """One-line JSON success receipt printed to stdout (no biological data)."""

    destination: str
    sha256: str


def _parse_sha256(raw: str) -> str:
    """Lowercased digest if it is exactly 64 hex characters, else UsageError."""
    if not _SHA256_SHAPE.fullmatch(raw):
        raise UsageError(
            "--sha256 must be exactly 64 hex characters",
            code="USAGE_ERROR",
            context={"sha256": raw},
        )
    return raw.lower()


def install_verified(source: Path, target: Path, expected: str) -> FetchReceipt:
    """Copy source to target with a streaming SHA256 check; atomic on success.

    Bytes land in a NamedTemporaryFile in the target's parent (same
    filesystem) and are os.replace()d over the target only after the digest
    verifies, so any failure leaves an existing target untouched. The temp
    file is removed on every error path.
    """
    if not source.is_file():
        raise InputError(
            f"source file not found or unreadable: {source}",
            code="INPUT_NOT_FOUND",
            context={"source": str(source)},
        )
    temp_path: Path | None = None
    try:
        with source.open("rb") as src, NamedTemporaryFile(dir=target.parent, delete=False) as temp:
            temp_path = Path(temp.name)
            hasher = hashlib.sha256()
            while chunk := src.read(_CHUNK_BYTES):
                hasher.update(chunk)
                temp.write(chunk)
            actual = hasher.hexdigest()
            if actual != expected:
                raise InputError(
                    f"SHA256 mismatch for {source}",
                    code="CHECKSUM_MISMATCH",
                    context={
                        "source": str(source),
                        "expected": expected,
                        "actual": actual,
                    },
                )
        assert temp_path is not None
        os.replace(temp_path, target)
    except OSError as exc:
        raise InputError(
            f"cannot install {source} to {target}: {exc}",
            code="FILE_IO_ERROR",
            context={"source": str(source), "target": str(target)},
        ) from exc
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
    return FetchReceipt(destination=str(target), sha256=expected)


def db_install_command(
    source: Annotated[
        Path,
        typer.Argument(
            help="Local source file path (plain filesystem only; no URLs, no provider IDs)."
        ),
    ],
    sha256: Annotated[
        str,
        typer.Option("--sha256", help="Expected SHA256 digest of the source (64 hex chars)."),
    ],
    output: Annotated[
        Path,
        typer.Option(
            "--output", help="Destination path; replaced atomically only after verification."
        ),
    ],
) -> None:
    """Install a local file to --output after verifying its SHA256.

    VERIFIED LOCAL-FILE INSTALLATION ONLY: this never fetches over the
    network and knows nothing about database providers or sequence content —
    it checksum-verifies and atomically installs bytes. On success a one-line
    JSON receipt (destination, verified digest) is printed to stdout.
    """

    def run() -> None:
        expected = _parse_sha256(sha256)
        receipt = install_verified(source, output, expected)
        typer.echo(receipt.model_dump_json())

    dispatch(run)
