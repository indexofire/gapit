"""Streaming FASTA reader (plain, gzip, bzip2)."""

import bz2
import gzip
from collections.abc import Iterator
from pathlib import Path
from typing import IO

from pydantic import BaseModel

from gapit.errors import InputError


class FastaRecord(BaseModel, frozen=True):
    """One FASTA record: ``id`` is the first whitespace-delimited header token."""

    id: str
    description: str
    sequence: str


def _open_text(path: Path) -> IO[str]:
    """Open a FASTA file as UTF-8 text, transparently decompressing .gz/.bz2."""
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    if path.name.endswith(".bz2"):
        return bz2.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def _finalize(path: Path, pending: tuple[str, str, list[str]]) -> FastaRecord:
    """Turn accumulated (id, description, sequence lines) into a FastaRecord."""
    seqid, description, chunks = pending
    if not chunks:
        raise InputError(
            f"{path}: record {seqid!r} has an empty sequence",
            code="INVALID_FASTA",
            context={"file": str(path)},
        )
    return FastaRecord(id=seqid, description=description, sequence="".join(chunks))


def iter_fasta(path: Path) -> Iterator[FastaRecord]:
    """Stream FastaRecords from a FASTA file (plain/.gz/.bz2, UTF-8).

    Raises InputError on content before the first ``>`` header or on a record
    with an empty sequence; an empty file yields zero records.
    """
    with _open_text(path) as handle:
        # pending = (id, description, sequence lines) of the record being read
        pending: tuple[str, str, list[str]] | None = None
        for lineno, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if not line.startswith(">"):
                if pending is None:
                    raise InputError(
                        f"{path}: content before first '>' header at line {lineno}",
                        code="INVALID_FASTA",
                        context={"file": str(path)},
                    )
                pending[2].append(line)
                continue
            if pending is not None:
                yield _finalize(path, pending)
            parts = line[1:].strip().split(maxsplit=1)
            pending = (
                parts[0] if parts else "",
                parts[1] if len(parts) > 1 else "",
                [],
            )
        if pending is not None:
            yield _finalize(path, pending)
