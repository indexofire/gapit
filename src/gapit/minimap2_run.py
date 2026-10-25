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
from collections.abc import Iterator
from pathlib import Path

from gapit.db import Database
from gapit.errors import DependencyError, GapitError
from gapit.gseek_engine import engine, iter_paf_gseek
from gapit.paf import PafRecord, Preset, parse_paf_row


def _iter_minimap2(argv: list[str], r1: Path) -> Iterator[PafRecord]:
    """Run one minimap2 invocation, yielding PAF rows as stdout lines arrive
    (no full-PAF materialization). stderr is drained by a background thread:
    an undrained stderr pipe fills (~64 KB) and blocks the child mid-run.
    Returncode != 0 raises MINIMAP2_FAILED with the captured stderr text."""
    try:
        process = subprocess.Popen(
            argv,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
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
    try:
        if process.stdout is not None:
            for line in process.stdout:
                stripped = line.rstrip("\n")
                if stripped.strip():
                    yield parse_paf_row(stripped)
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


def _iter_legacy(
    lanes: list[tuple[Path, Path | None]],
    database: Database,
    *,
    read_type: Preset,
    threads: int,
    debug: bool = False,
    nm_tags: bool = False,
) -> Iterator[PafRecord]:
    """The original subprocess path; kept as the legacy/shadow oracle."""
    for r1, r2 in lanes:
        argv = [
            "minimap2",
            "-x",
            read_type,
            "-t",
            str(threads),
        ]
        if read_type == "sr":
            # Diverged-allele sensitivity for short reads: the stock sr seed
            # (k=21) finds no exact run on a ~90%-identity 117 bp marker gene
            # (26ECO0071 astA: 0 alignments at k=21, 55% primary breadth at
            # k=13, 100% at k=11), silently flipping typing calls. k=11 + a
            # softer mismatch cost recover it; presence stays breadth-gated.
            argv += ["-k11", "-w6", "-B2"]
        if nm_tags:
            argv.append("--cs")
        # Input paths go to argv absolutized: an absolute path starts with
        # "/" and can never parse as a minimap2 option, so a query file
        # named "-d" cannot make minimap2 dump an index over its mate path.
        argv += [str(database.sequences_path.absolute()), str(r1.absolute())]
        if r2 is not None:
            argv.append(str(r2.absolute()))
        if debug:
            print(f"gapit: run: {shlex.join(argv)}", file=sys.stderr)
        yield from _iter_minimap2(argv, r1)


def iter_minimap2(
    lanes: list[tuple[Path, Path | None]],
    database: Database,
    *,
    read_type: Preset,
    threads: int,
    debug: bool = False,
    nm_tags: bool = False,
) -> Iterator[PafRecord]:
    """Stream PAF rows from one minimap2 invocation per lane.

    The public iterator is the production path (reads/cluster consume it
    incrementally); ``run_minimap2`` remains the list-returning convenience
    wrapper for callers/tests that want a materialized view. In gseek mode the
    result is materialized per lane by gseek's alignment API.
    """
    switch = engine()
    if switch.mode.value == "legacy":
        yield from _iter_legacy(
            lanes,
            database,
            read_type=read_type,
            threads=threads,
            debug=debug,
            nm_tags=nm_tags,
        )
        return
    if switch.mode.value == "gseek":
        yield from iter_paf_gseek(
            lanes,
            database,
            read_type=read_type,
            threads=threads,
            debug=debug,
            nm_tags=nm_tags,
        )
        return
    yield from switch.run(
        lambda: list(
            _iter_legacy(
                lanes,
                database,
                read_type=read_type,
                threads=threads,
                debug=debug,
                nm_tags=nm_tags,
            )
        ),
        lambda: list(
            iter_paf_gseek(
                lanes,
                database,
                read_type=read_type,
                threads=threads,
                debug=debug,
                nm_tags=nm_tags,
            )
        ),
    )


def run_minimap2(
    lanes: list[tuple[Path, Path | None]],
    database: Database,
    *,
    read_type: Preset,
    threads: int,
    debug: bool = False,
    nm_tags: bool = False,
) -> list[PafRecord]:
    """Materialized convenience view over :func:`iter_minimap2`; production
    reads/cluster paths consume the iterator to keep memory O(targets)."""
    return list(
        iter_minimap2(
            lanes,
            database,
            read_type=read_type,
            threads=threads,
            debug=debug,
            nm_tags=nm_tags,
        )
    )
