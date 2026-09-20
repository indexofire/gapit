"""The minimap2 invocation layer for reads mode (SPEC.md §10).

Split from reads.py (which was over the 250 pure-LOC ceiling) so this module
owns exactly one thing: running minimap2 and streaming its PAF rows. The
minimap2 invocation's own preset vocabulary (``ReadType``) lives in paf.py,
the lowest layer of the reads stack. Not to be confused with minimap.py,
which owns the blastn-side COVERAGE_MAP arithmetic.
"""

import shlex
import subprocess
import sys
import threading
from pathlib import Path

from gapit.db import Database
from gapit.errors import DependencyError, GapitError
from gapit.paf import PafRecord, ReadType, parse_paf_row


def _stream_minimap2(argv: list[str], r1: Path) -> list[PafRecord]:
    """Run one minimap2 invocation, streaming PAF rows off stdout as lines
    arrive (no full-PAF materialization). stderr is drained by a background
    thread: an undrained stderr pipe fills (~64 KB) and blocks the child
    mid-run. Returncode != 0 raises MINIMAP2_FAILED with the captured
    stderr text."""
    try:
        process = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except FileNotFoundError as exc:
        raise DependencyError(
            "required binary not found on PATH: minimap2",
            code="MISSING_DEPENDENCY",
            context={"binary": "minimap2"},
        ) from exc
    stderr_chunks: list[str] = []

    def drain() -> None:
        if process.stderr is not None:
            stderr_chunks.append(process.stderr.read())

    thread = threading.Thread(target=drain)
    thread.start()
    rows: list[PafRecord] = []
    try:
        if process.stdout is not None:
            for line in process.stdout:
                stripped = line.rstrip("\n")
                if stripped.strip():
                    rows.append(parse_paf_row(stripped))
    finally:
        if process.stdout is not None:
            process.stdout.close()  # on a parse error, SIGPIPE stops the child
        process.wait()
        thread.join()
    if process.returncode != 0:
        raise GapitError(
            f"minimap2 failed: {''.join(stderr_chunks).strip()}",
            code="MINIMAP2_FAILED",
            context={"binary": "minimap2", "file": str(r1)},
        )
    return rows


def run_minimap2(
    lanes: list[tuple[Path, Path | None]],
    database: Database,
    *,
    read_type: ReadType,
    threads: int,
    debug: bool = False,
    nm_tags: bool = False,
) -> list[PafRecord]:
    """Run one minimap2 invocation per lane (PAF on stdout) and concatenate
    the rows, streaming each lane's output row by row. The index argument is
    always the ``sequences`` FASTA — minimap2 indexes it in memory with the
    invocation preset's own parameters. A persisted ``.mmi`` is deliberately
    rejected even when one sits beside the FASTA: a default-built index
    overrides the ``-x`` preset's indexing parameters (``-k, -w or -H
    overridden by prebuilt index``), which misassigns close homologs and
    benchmarks slower than in-memory indexing (2026-09-19: blaCTX-M/blaSHV
    allele divergence, +1.2 s on the ncbi db). minimap2's pairing semantics
    for >2 input files are undocumented; per-lane runs (r1[i] alone or with
    its mate r2[i]) are deterministic. With ``nm_tags``, ``--cs`` is added so
    rows carry ``NM:i:`` (minimap2 omits NM from PAF output without it; the
    alignments themselves are unchanged) — used by gapit.reads/2 identity
    filtering. With ``debug``, echo each argv to stderr (abricate --debug
    parity)."""
    rows: list[PafRecord] = []
    for r1, r2 in lanes:
        argv = [
            "minimap2",
            "-x",
            read_type,
            "-t",
            str(threads),
        ]
        if nm_tags:
            argv.append("--cs")
        argv += [str(database.sequences_path), str(r1)]
        if r2 is not None:
            argv.append(str(r2))
        if debug:
            print(f"gapit: run: {shlex.join(argv)}", file=sys.stderr)
        rows.extend(_stream_minimap2(argv, r1))
    return rows
