"""Batched gene screening through gseek.align_many.

The batch is an optimization layer over the normal per-file gene path: only
all-FASTA inputs with a gseek build exposing ``annotate_batch`` are eligible.
Records are split by gseek's ``batch_index`` tag, then each file goes through
the same ``report_from_rows`` fold as the sequential path. With ``--jobs`` the
inputs are split into at most that many contiguous chunks; each chunk is one
aligner invocation, and chunk results are reassembled in input order.
"""

import math
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Literal

from gapit import db
from gapit.blast import report_from_rows
from gapit.errors import InputError
from gapit.gseek_engine import batch_supported, engine_mode, run_blast_many_gseek
from gapit.report import Report, ScreeningParams
from gapit.seqconvert import SeqFormat, detect_format


def _all_fasta(paths: Sequence[Path]) -> bool:
    try:
        return all(detect_format(path) is SeqFormat.fasta for path in paths)
    except InputError:
        # Let the per-file path produce the documented invalid-input envelope.
        return False


def screen_batch(
    paths: list[Path],
    database: db.Database,
    params: ScreeningParams,
    *,
    dbtype: Literal["nucl", "prot"],
    debug: bool,
    merge_fragments: bool,
    jobs: int,
) -> list[Report] | None:
    """Return reports for a gseek batch, or None when batching is unavailable."""
    if engine_mode() != "gseek" or not batch_supported() or not _all_fasta(paths):
        return None
    chunk_size = max(1, math.ceil(len(paths) / max(1, jobs)))
    chunks = [paths[start : start + chunk_size] for start in range(0, len(paths), chunk_size)]

    def run_chunk(chunk: list[Path]) -> list[Report] | None:
        groups = run_blast_many_gseek(
            chunk,
            database,
            params,
            dbtype=dbtype,
            debug=debug,
        )
        if groups is None:
            return None
        return [
            report_from_rows(rows, path, database, params, merge_fragments=merge_fragments)
            for path, rows in zip(chunk, groups, strict=True)
        ]

    if jobs == 1 or len(chunks) == 1:
        results = [run_chunk(chunk) for chunk in chunks]
    else:
        with ThreadPoolExecutor(max_workers=min(jobs, len(chunks))) as executor:
            results = list(executor.map(run_chunk, chunks))
    if any(result is None for result in results):
        return None
    return [report for result in results for report in result or ()]


__all__ = ["screen_batch"]
