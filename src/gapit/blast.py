"""BLAST invocation, tabular parsing, and the screening pipeline (SPEC.md §3)."""

import re
import shlex
import subprocess
import sys
import threading
from pathlib import Path
from typing import Literal

from pydantic import BaseModel

from gapit.db import Database
from gapit.errors import DependencyError, GapitError, InputError
from gapit.hits import process_rows
from gapit.report import Report, ScreeningParams

BLAST_FIELDS = [
    "qseqid",
    "qstart",
    "qend",
    "qlen",
    "sseqid",
    "sstart",
    "send",
    "slen",
    "sstrand",
    "evalue",
    "length",
    "pident",
    "gaps",
    "gapopen",
    "stitle",
]

_OUTFMT = "6 " + " ".join(BLAST_FIELDS)
_MIN_BLAST_VERSION = (2, 2, 30)


class BlastRow(BaseModel, frozen=True):
    """One outfmt-6 row, typed at the boundary (SPEC.md §3 field order)."""

    qseqid: str
    qstart: int
    qend: int
    qlen: int
    sseqid: str
    sstart: int
    send: int
    slen: int
    sstrand: str
    evalue: float
    length: int
    pident: float
    gaps: int
    gapopen: int
    stitle: str


def parse_blast_row(line: str) -> BlastRow:
    """Parse one tab-delimited outfmt-6 line; a row with != 15 columns is a
    hard error (upstream wording)."""
    fields = line.split("\t")
    if len(fields) != 15:
        raise GapitError("can not find sequence data", code="BLAST_PARSE_FAILED")
    return BlastRow(
        qseqid=fields[0],
        qstart=int(fields[1]),
        qend=int(fields[2]),
        qlen=int(fields[3]),
        sseqid=fields[4],
        sstart=int(fields[5]),
        send=int(fields[6]),
        slen=int(fields[7]),
        sstrand=fields[8],
        evalue=float(fields[9]),
        length=int(fields[10]),
        pident=float(fields[11]),
        gaps=int(fields[12]),
        gapopen=int(fields[13]),
        stitle=fields[14],
    )


def ensure_blast() -> None:
    """Require blastn >= 2.2.30 on PATH (SPEC.md §1 version gate)."""
    try:
        result = subprocess.run(["blastn", "-version"], check=False, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise DependencyError(
            "required binary not found on PATH: blastn",
            code="MISSING_DEPENDENCY",
            context={"binary": "blastn"},
        ) from exc
    if result.returncode != 0:
        raise DependencyError(
            f"blastn -version failed: {result.stderr.strip()}",
            code="BLAST_VERSION_CHECK_FAILED",
        )
    match = re.search(r"blastn:\s*(\d+)\.(\d+)\.(\d+)", result.stdout)
    if match is None:
        raise DependencyError(
            f"could not parse blastn version from: {result.stdout!r}",
            code="BLAST_VERSION_CHECK_FAILED",
        )
    major, minor, revision = (int(part) for part in match.groups())
    if (major, minor, revision) < _MIN_BLAST_VERSION:
        raise DependencyError(
            f"blastn {major}.{minor}.{revision} is older than 2.2.30",
            code="BLAST_VERSION_TOO_OLD",
        )


def _pipeline(query: Path, argv: list[str]) -> str:
    """`any2fasta -q -u <query> | <argv>`; returns blast's stdout as text."""
    try:
        any2fasta = subprocess.Popen(
            ["any2fasta", "-q", "-u", str(query)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError as exc:
        raise DependencyError(
            "required binary not found on PATH: any2fasta",
            code="MISSING_DEPENDENCY",
            context={"binary": "any2fasta"},
        ) from exc
    try:
        blast = subprocess.Popen(
            argv, stdin=any2fasta.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
    except FileNotFoundError as exc:
        raise DependencyError(
            f"required binary not found on PATH: {argv[0]}",
            code="MISSING_DEPENDENCY",
            context={"binary": argv[0]},
        ) from exc
    if any2fasta.stdout is not None:
        any2fasta.stdout.close()  # blast owns the read end now; SIGPIPE propagates
    # stderr is drained by a background thread: an undrained stderr pipe
    # fills (~64 KB) and blocks any2fasta mid-run (minimap2_run._stream_minimap2).
    any2fasta_err_chunks: list[bytes] = []

    def drain() -> None:
        if any2fasta.stderr is not None:
            any2fasta_err_chunks.append(any2fasta.stderr.read())

    drain_thread = threading.Thread(target=drain)
    drain_thread.start()
    try:
        blast_out, blast_err = blast.communicate()
        any2fasta_rc = any2fasta.wait()
    finally:
        drain_thread.join()
    any2fasta_err = b"".join(any2fasta_err_chunks)
    # blast first: if it crashed, its stderr names the real cause (any2fasta
    # may merely have taken the SIGPIPE).
    if blast.returncode != 0:
        raise GapitError(
            f"{argv[0]} failed: {blast_err.decode('utf-8', 'replace').strip()}",
            code="BLAST_FAILED",
            context={"binary": argv[0]},
        )
    if any2fasta_rc != 0:
        raise InputError(
            f"invalid input file {query}: {any2fasta_err.decode('utf-8', 'replace').strip()}",
            code="INVALID_INPUT",
            context={"file": str(query)},
        )
    return blast_out.decode("utf-8")


def run_screen(
    query: Path,
    database: Database,
    params: ScreeningParams,
    *,
    dbtype: Literal["nucl", "prot"],
    debug: bool = False,
) -> list[BlastRow]:
    """Run the any2fasta -> blastn/blastx pipeline for one query file.

    The database must already be indexed; protein databases switch to blastx
    without -perc_identity (upstream quirk, minid silently ignored). With
    ``debug``, echo the exact any2fasta and blast argv to stderr
    (abricate --debug parity). ``dbtype`` comes from the caller resolving it
    once per run, so dependency/index errors surface before any per-file work.
    """
    if dbtype == "prot":
        argv = [
            "blastx",
            "-task",
            "blastx-fast",
            "-seg",
            "no",
            "-db",
            str(database.sequences_path),
            "-outfmt",
            _OUTFMT,
            "-num_threads",
            str(params.threads),
            "-evalue",
            "1E-20",
            "-culling_limit",
            "1",
            "-max_target_seqs",
            "10000",
        ]
        sys.stderr.write("--minid is not applied to protein databases (abricate parity)\n")
    else:
        argv = [
            "blastn",
            "-task",
            "blastn",
            "-dust",
            "no",
            "-perc_identity",
            str(params.minid),
            "-db",
            str(database.sequences_path),
            "-outfmt",
            _OUTFMT,
            "-num_threads",
            str(params.threads),
            "-evalue",
            "1E-20",
            "-culling_limit",
            "1",
            "-max_target_seqs",
            "10000",
        ]
    if debug:
        sys.stderr.write(f"gapit: run: {shlex.join(['any2fasta', '-q', '-u', str(query)])}\n")
        sys.stderr.write(f"gapit: run: {shlex.join(argv)}\n")
    output = _pipeline(query, argv)
    return [parse_blast_row(line) for line in output.splitlines() if line.strip()]


def screen_file(
    query: Path,
    database: Database,
    params: ScreeningParams,
    *,
    dbtype: Literal["nucl", "prot"],
    debug: bool = False,
) -> Report:
    """Screen one input file against one database into a sorted Report.

    The caller owns the per-run gates (``ensure_blast`` and the one-shot
    ``dbtype`` resolution) so a multi-file run pays each probe once."""
    rows = run_screen(query, database, params, dbtype=dbtype, debug=debug)
    hits = process_rows(rows, mincov=params.mincov, default_db=params.db)
    ordered = sorted(hits, key=lambda hit: (hit.sequence, hit.start))
    return Report(file=str(query), hits=tuple(ordered))
