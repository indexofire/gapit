"""Native input normalization: FASTA/FASTQ/GenBank/EMBL (.gz/.bz2) → FASTA lines.

Replaces the external ``any2fasta -q -u <file>`` stage of the screening
pipeline (blast.py). Semantics extracted from any2fasta 0.8.1 (Perl,
bioconda; ``.pixi/envs/default/bin/any2fasta``) invoked with ``-q -u`` only —
no ``-n`` (N-purification), ``-l``, ``-g`` (VERSION), or ``-s`` (description
stripping). Perl line references below are that file.

- Detection (L89-99, L136-153): first line of the decompressed stream, in
  order GENBANK ``^LOCUS\\h``, EMBL ``^ID\\h``, FASTA ``^>\\S``, FASTQ
  ``^@\\S`` (``\\h`` = space/tab; ``\\S`` = non-whitespace after the marker).
  Empty input dies with "The input appears to be empty" (L131-133); an
  unrecognized first line dies with "Unfamilar format with first line: ..."
  (L153 — typo theirs, kept for message parity with the retired binary).
- purify_dna (L164-170): with ``-u`` only uc() — sequences uppercased, nothing
  else touched. purify_id (L174-181) is a no-op without ``-s``: headers pass
  through verbatim.
- FASTA (L185-200): header lines printed verbatim; sequence lines uppercased;
  blank/whitespace-only lines skipped entirely (L189).
- FASTQ (L204-215): strict 4-line stride. Header = the ``@`` line minus its
  leading ``@``, verbatim (id + description); sequence = uc() of line 2;
  ``+``/quality lines dropped. A trailing lone ``@`` header emits nothing
  (loop guard ``$i < $#lines``, L208).
- GenBank (L255-296): id = first token after ``LOCUS\\s+`` (L285-287).
  DEFINITION is never read — records carry no description. VERSION overrides
  the id only under ``-g`` (L290-292; unused). ORIGIN data lines: drop the
  first 10 columns (coordinate prefix, L279), then strip whitespace (L280) —
  digits past column 10 are kept. Records flush at ``//`` (L264-269); a
  record never terminated by ``//`` is dropped.
- EMBL (L300-339): id = text after ``ID\\s+`` up to the first ``;`` (L330-331,
  not trimmed). DE is never read. SQ data lines: strip all whitespace AND
  digits (L325). Flush at ``//``.

Deliberate divergences, invisible to BLAST: sequence is re-wrapped at 60
columns (the Perl reuses input wrapping) and input is read with universal
newlines (the Perl preserves ``\\r``). Parsed records — id, description,
uppercased sequence — are identical, so the blastn query is unchanged and
abricate parity is unaffected. gapit owns no any2fasta dependency; the
differential suite runs the real binary from the parity env's PATH, where
abricate provides it transitively (tests/test_seqconvert_differential.py).
"""

import re
from collections.abc import Iterator
from enum import Enum
from pathlib import Path
from typing import IO

from gapit.errors import InputError
from gapit.fasta import open_text

_WRAP = 60


class SeqFormat(Enum):
    """Input formats detected from the first decompressed line."""

    fasta = "fasta"
    fastq = "fastq"
    genbank = "genbank"
    embl = "embl"


def _tagged(line: str, tag: str) -> bool:
    """``^tag\\h`` — tag at column 1 followed by one space or tab."""
    return line.startswith(tag) and line[len(tag) : len(tag) + 1] in (" ", "\t")


def detect_format(path: Path) -> SeqFormat:
    """Content-sniff the first decompressed line (reads.py detect_read_kind style)."""
    try:
        with open_text(path) as handle:
            first = handle.readline()
    except (OSError, UnicodeDecodeError, EOFError) as exc:
        raise InputError(
            f"could not read input: {exc}",
            code="INVALID_INPUT",
            context={"file": str(path)},
        ) from exc
    if not first:
        raise InputError(
            "The input appears to be empty",
            code="INVALID_INPUT",
            context={"file": str(path)},
        )
    line = first.rstrip("\r\n")
    if _tagged(line, "LOCUS"):
        return SeqFormat.genbank
    if _tagged(line, "ID"):
        return SeqFormat.embl
    if len(line) > 1 and line[0] == ">" and not line[1].isspace():
        return SeqFormat.fasta
    if len(line) > 1 and line[0] == "@" and not line[1].isspace():
        return SeqFormat.fastq
    raise InputError(
        f"Unfamilar format with first line: {line}",
        code="INVALID_INPUT",
        context={"file": str(path)},
    )


def _iter_fasta(handle: IO[str]) -> Iterator[tuple[str, str]]:
    """(header, sequence) pairs; headers verbatim, blank lines skipped."""
    header: str | None = None
    chunks: list[str] = []
    for raw in handle:
        line = raw.rstrip("\n")
        if not line.strip():
            continue
        if line.startswith(">"):
            if header is not None:
                yield (header, "".join(chunks).upper())
            header, chunks = line[1:], []
        else:
            chunks.append(line)
    if header is not None:
        yield (header, "".join(chunks).upper())


def _iter_fastq(handle: IO[str]) -> Iterator[tuple[str, str]]:
    """Stride-4 records; a trailing lone '@' header emits nothing (perl guard)."""
    lines = iter(handle)
    for header_line in lines:
        sequence_line = next(lines, None)
        if sequence_line is None:
            return
        yield (header_line.rstrip("\n")[1:], sequence_line.rstrip("\n").upper())
        next(lines, None)  # '+' separator
        next(lines, None)  # quality line


def _iter_genbank(handle: IO[str]) -> Iterator[tuple[str, str]]:
    """LOCUS-name ids without description; ORIGIN columns 11+ minus whitespace."""
    acc = ""
    chunks: list[str] = []
    in_seq = False
    for raw in handle:
        line = raw.rstrip("\n")
        if line.startswith("//"):
            yield (acc, "".join(chunks).upper())
            acc, chunks, in_seq = "", [], False
        elif line.startswith("ORIGIN"):
            in_seq = True
        elif in_seq:
            chunks.append(re.sub(r"\s+", "", line[10:]))
        elif (match := re.match(r"LOCUS\s+(\S+)", line)) is not None:
            acc = match.group(1)


def _iter_embl(handle: IO[str]) -> Iterator[tuple[str, str]]:
    """ID-before-';' ids without description; SQ lines minus whitespace+digits."""
    acc = ""
    chunks: list[str] = []
    in_seq = False
    for raw in handle:
        line = raw.rstrip("\n")
        if line.startswith("//"):
            yield (acc, "".join(chunks).upper())
            acc, chunks, in_seq = "", [], False
        elif re.match(r"SQ\s", line) is not None:
            in_seq = True
        elif in_seq:
            chunks.append(re.sub(r"[\s\d]", "", line))
        elif (match := re.match(r"ID\s+([^;]+)", line)) is not None:
            acc = match.group(1)


def _records(fmt: SeqFormat, handle: IO[str]) -> Iterator[tuple[str, str]]:
    records: Iterator[tuple[str, str]]
    match fmt:
        case SeqFormat.fasta:
            records = _iter_fasta(handle)
        case SeqFormat.fastq:
            records = _iter_fastq(handle)
        case SeqFormat.genbank:
            records = _iter_genbank(handle)
        case SeqFormat.embl:
            records = _iter_embl(handle)
    return records


def to_fasta_lines(path: Path, fmt: SeqFormat | None = None) -> Iterator[str]:
    """Yield complete FASTA lines for one input file: one verbatim header line
    per record, then its uppercased sequence wrapped at 60 columns.

    ``fmt`` may be pre-supplied to reuse a sniff (blast.run_blastn's debug
    echo shares one detection). Raises InputError (INVALID_INPUT) on unreadable,
    empty, or unrecognized input — the retired any2fasta's fatal path.
    """
    if fmt is None:
        fmt = detect_format(path)
    try:
        with open_text(path) as handle:
            for header, sequence in _records(fmt, handle):
                yield f">{header}\n"
                for offset in range(0, len(sequence), _WRAP):
                    yield f"{sequence[offset : offset + _WRAP]}\n"
    except (OSError, UnicodeDecodeError, EOFError) as exc:
        raise InputError(
            f"could not read input: {exc}",
            code="INVALID_INPUT",
            context={"file": str(path)},
        ) from exc
