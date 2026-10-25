"""Optional gseek bridge for BLAST and minimap2 execution.

``GAPIT_ENGINE`` selects the implementation:

- ``legacy`` (default): the existing subprocess paths, no gseek import;
- ``gseek``: gseek plans and parses the same BLAST/minimap2 work;
- ``shadow``: run both, log coarse differences, and return the legacy result.

gseek is optional and requires Python >= 3.12. Importing this module never
imports gseek; only a ``gseek``/``shadow`` call does.
"""

from __future__ import annotations

import importlib
import logging
import os
import sys
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Literal, Protocol, cast

from gapit.db import Database
from gapit.errors import DependencyError, GapitError
from gapit.report import ScreeningParams

if TYPE_CHECKING:
    from gapit.blast import BlastRow
    from gapit.paf import PafRecord

_LOGGER = logging.getLogger("gapit.gseek")
_ENGINE_ENV = "GAPIT_ENGINE"
_VALID_MODES = frozenset({"legacy", "gseek", "shadow"})


class _GseekRecord(Protocol):
    """The subset of a gseek AlignmentRecord used by gapit."""

    query_id: str
    target_id: str
    query_length: int | None
    target_length: int | None
    query_start: int | None
    query_end: int | None
    target_start: int | None
    target_end: int | None
    identity: float | None
    aligned_length: int | None
    matches: int | None
    gap_count: int | None
    gap_openings: int | None
    evalue: float | None
    mapq: int | None
    target_description: str | None
    is_primary: bool
    tags: Mapping[str, str]

    @property
    def is_reverse(self) -> bool: ...


class _GseekResult(Protocol):
    records: tuple[_GseekRecord, ...]


class _GseekModule(Protocol):
    def align(self, **kwargs: object) -> _GseekResult: ...


def engine_mode() -> str:
    """Return ``GAPIT_ENGINE`` (``legacy`` by default), validating its value."""
    raw = os.environ.get(_ENGINE_ENV, "legacy").strip().lower() or "legacy"
    if raw not in _VALID_MODES:
        raise GapitError(
            f"unknown {_ENGINE_ENV}={raw!r}; expected one of {sorted(_VALID_MODES)}",
            code="INVALID_ENGINE",
            context={_ENGINE_ENV: raw},
        )
    return raw


def _gseek() -> _GseekModule:
    try:
        module = importlib.import_module("gseek")
    except ImportError as exc:
        raise DependencyError(
            f"{_ENGINE_ENV}=gseek requires the optional gseek package "
            "(Python >= 3.12; install gapit[gseek])",
            code="MISSING_DEPENDENCY",
            context={"package": "gseek"},
        ) from exc
    return cast(_GseekModule, module)


def _record_int(value: int | None, field: str) -> int:
    if value is None:
        raise GapitError(
            f"gseek {field} field is missing from an alignment record",
            code="GSEEK_RECORD_INVALID",
            context={"field": field},
        )
    return value


def _blast_row(record: _GseekRecord, *, protein: bool = False) -> BlastRow:
    """Map one gseek BLAST record back to BLAST's 15-column row contract.

    gseek stores 0-based half-open, forward-strand coordinates; BLAST emits
    1-based inclusive coordinates and reports a minus subject as sstart > send.
    gapit's dedup/swap logic depends on exactly that legacy representation.
    """
    from gapit.blast import BlastRow

    qstart = _record_int(record.query_start, "query_start") + 1
    qend = _record_int(record.query_end, "query_end")
    tstart = _record_int(record.target_start, "target_start")
    tend = _record_int(record.target_end, "target_end")
    reverse = record.is_reverse
    sstart = tend if reverse else tstart + 1
    send = tstart + 1 if reverse else tend
    return BlastRow(
        qseqid=record.query_id,
        qstart=qstart,
        qend=qend,
        qlen=_record_int(record.query_length, "query_length"),
        sseqid=record.target_id,
        sstart=sstart,
        send=send,
        slen=_record_int(record.target_length, "target_length"),
        sstrand="N/A" if protein else ("minus" if reverse else "plus"),
        evalue=record.evalue if record.evalue is not None else 0.0,
        length=_record_int(record.aligned_length, "aligned_length"),
        pident=record.identity if record.identity is not None else 0.0,
        gaps=record.gap_count or 0,
        gapopen=record.gap_openings or 0,
        stitle=record.target_description or record.target_id,
    )


_INDEX_POLICY: dict[Path, str] = {}
_MINIMAP2_NAME_LIMIT = 255


def _minimap2_index_policy(database: Database) -> str:
    """Use a direct FASTA target when minimap2 names exceed its index limit.

    minimap2's ``.mmi`` index replaces over-long target names with numeric
    sequence indices; the legacy path passes the FASTA target directly and
    keeps the names. CARD has six such names (including MexB), so an index
    changes the reported gene even though the alignment is identical.
    """
    path = database.sequences_path
    cached = _INDEX_POLICY.get(path)
    if cached is not None:
        return cached
    policy = "auto"
    if not path.is_file():
        return policy
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split(maxsplit=1)[0]
                if len(name) > _MINIMAP2_NAME_LIMIT:
                    policy = "none"
                    break
    _INDEX_POLICY[path] = policy
    return policy


def _paf_record(record: _GseekRecord) -> PafRecord:
    """Map one gseek minimap2 record to gapit's PAF boundary model."""
    from gapit.paf import PafRecord

    tags = record.tags
    nm_raw = tags.get("NM")
    try:
        nm = int(nm_raw) if nm_raw is not None else None
    except ValueError as exc:
        raise GapitError(
            f"gseek minimap2 NM tag is not an integer: {nm_raw!r}",
            code="GSEEK_RECORD_INVALID",
            context={"tag": "NM"},
        ) from exc
    return PafRecord(
        qname=record.query_id,
        qlen=_record_int(record.query_length, "query_length"),
        qstart=_record_int(record.query_start, "query_start"),
        qend=_record_int(record.query_end, "query_end"),
        strand="-" if record.is_reverse else "+",
        tname=record.target_id,
        tlen=_record_int(record.target_length, "target_length"),
        tstart=_record_int(record.target_start, "target_start"),
        tend=_record_int(record.target_end, "target_end"),
        nmatch=_record_int(record.matches, "matches"),
        alen=_record_int(record.aligned_length, "aligned_length"),
        mapq=_record_int(record.mapq, "mapq"),
        is_primary=record.is_primary,
        nm=nm,
        cs=tags.get("cs"),
    )


def run_blast_gseek(
    fasta_text: str,
    database: Database,
    params: ScreeningParams,
    *,
    dbtype: Literal["nucl", "prot"],
    debug: bool = False,
) -> list[BlastRow]:
    """Run one normalized query through gseek's BLAST adapter."""
    max_targets = 10000
    gseek_params: dict[str, object] = {
        "task": "blastx-fast" if dbtype == "prot" else "blastn",
        "evalue": 1e-20,
        "max_target_seqs": max_targets,
        "culling_limit": 1,
    }
    if dbtype == "prot":
        gseek_params["seg"] = False
    else:
        gseek_params["dust"] = False
        gseek_params["min_identity"] = params.minid
    if debug:
        print(
            f"gapit: gseek blast: db={database.name} dbtype={dbtype}",
            file=sys.stderr,
        )
    with tempfile.TemporaryDirectory(prefix="gapit-gseek-") as tmp:
        query_path = Path(tmp) / "query.fasta"
        query_path.write_text(fasta_text, encoding="utf-8")
        gseek = _gseek()
        try:
            result = gseek.align(
                backend="blast",
                query=query_path,
                target=database.sequences_path,
                query_kind="fasta",
                params=gseek_params,
                threads=params.threads,
                index="auto",
                record_manifest=False,
            )
        except Exception as exc:
            raise GapitError(
                f"gseek BLAST failed: {exc}",
                code="GSEEK_FAILED",
                context={"backend": "blast"},
            ) from exc
    return [_blast_row(record, protein=dbtype == "prot") for record in result.records]


def iter_paf_gseek(
    lanes: Sequence[tuple[Path, Path | None]],
    database: Database,
    *,
    read_type: str,
    threads: int,
    debug: bool = False,
    nm_tags: bool = False,
) -> Iterable[PafRecord]:
    """Yield gapit PAF records from gseek minimap2 for every lane.

    ``secondary="all"`` means "do not pass --secondary=no": gapit's legacy
    invocation keeps minimap2's default secondary output and filters later.
    """
    gseek = _gseek()
    for r1, r2 in lanes:
        query = (r1,) if r2 is None else (r1, r2)
        gseek_params: dict[str, object] = {
            "preset": read_type,
            "secondary": "all",
            "cs": nm_tags,
        }
        if read_type == "sr":
            gseek_params.update({"k": 11, "w": 6, "B": 2})
        if debug:
            print(
                f"gapit: gseek minimap2: preset={read_type} reads={r1}",
                file=sys.stderr,
            )
        try:
            result = gseek.align(
                backend="minimap2",
                query=query,
                target=database.sequences_path,
                paired=r2 is not None,
                params=gseek_params,
                threads=threads,
                index=_minimap2_index_policy(database),
                record_manifest=False,
            )
        except Exception as exc:
            raise GapitError(
                f"gseek minimap2 failed: {exc}",
                code="GSEEK_FAILED",
                context={"backend": "minimap2"},
            ) from exc
        for record in result.records:
            yield _paf_record(record)


def log_shadow_error(engine: str, exc: BaseException) -> None:
    """Log a shadow-mode gseek failure without interrupting the legacy path."""
    _LOGGER.warning("gseek shadow %s failed; returning legacy: %s", engine, exc)


def log_blast_shadow(legacy: Sequence[BlastRow], new: Sequence[BlastRow]) -> None:
    """Log a coarse legacy-vs-gseek BLAST difference."""
    legacy_keys = {
        (row.qseqid, row.qstart, row.qend, row.sseqid, round(row.pident, 2)) for row in legacy
    }
    new_keys = {(row.qseqid, row.qstart, row.qend, row.sseqid, round(row.pident, 2)) for row in new}
    if legacy_keys != new_keys:
        _LOGGER.warning(
            "gseek shadow BLAST differs: legacy=%d gseek=%d only_legacy=%d only_gseek=%d",
            len(legacy),
            len(new),
            len(legacy_keys - new_keys),
            len(new_keys - legacy_keys),
        )


def log_paf_shadow(legacy: Sequence[PafRecord], new: Sequence[PafRecord]) -> None:
    """Log a coarse legacy-vs-gseek minimap2 difference."""
    legacy_keys = {
        (row.qname, row.tname, row.tstart, row.tend, row.strand, row.is_primary, row.mapq)
        for row in legacy
    }
    new_keys = {
        (row.qname, row.tname, row.tstart, row.tend, row.strand, row.is_primary, row.mapq)
        for row in new
    }
    if legacy_keys != new_keys:
        _LOGGER.warning(
            "gseek shadow minimap2 differs: legacy=%d gseek=%d only_legacy=%d only_gseek=%d",
            len(legacy),
            len(new),
            len(legacy_keys - new_keys),
            len(new_keys - legacy_keys),
        )


__all__ = [
    "engine_mode",
    "iter_paf_gseek",
    "log_blast_shadow",
    "log_paf_shadow",
    "log_shadow_error",
    "run_blast_gseek",
]
