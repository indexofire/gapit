"""PAF row parsing and interval arithmetic (minimap2 output boundary).

Typed parsing of minimap2's PAF rows (12 required fields + tags), the
half-open interval math used to aggregate alignment spans into coverage, and
the per-alignment identity rule behind gapit.reads/2 filtering.
"""

from collections.abc import Iterable
from typing import Literal

from pydantic import BaseModel

from gapit.errors import GapitError

ReadType = Literal["sr", "map-ont", "map-hifi"]


class PafRecord(BaseModel, frozen=True):
    """One PAF alignment row: 12 required fields + primary flag (tp:A:P) and
    mismatch count (NM:i:, None when minimap2 emitted no NM tag — it only
    does with ``--cs``)."""

    qname: str
    qlen: int
    qstart: int
    qend: int
    strand: Literal["+", "-"]
    tname: str
    tlen: int
    tstart: int
    tend: int
    nmatch: int
    alen: int
    mapq: int
    is_primary: bool = True
    nm: int | None = None


def parse_paf_row(line: str) -> PafRecord:
    """Parse one tab-delimited PAF line; <12 fields is a hard error. Records
    without a tp tag count as primary; NM:i: is extracted when present."""
    fields = line.split("\t")
    if len(fields) < 12:
        raise GapitError(
            f"can not parse PAF row (expected >= 12 fields): {line!r}",
            code="PAF_PARSE_FAILED",
        )
    is_primary = True
    nm: int | None = None
    for tag in fields[12:]:
        if tag.startswith("tp:A:"):
            is_primary = tag == "tp:A:P"
        elif tag.startswith("NM:i:"):
            nm = int(tag[len("NM:i:") :])
    strand = fields[4]
    if strand not in ("+", "-"):
        raise GapitError(
            f"can not parse PAF strand (expected + or -): {strand!r}",
            code="PAF_PARSE_FAILED",
        )
    # model_construct: the manual int()/strand checks above already guarantee
    # the field types; re-validating per PAF row would be redundant.
    return PafRecord.model_construct(
        qname=fields[0],
        qlen=int(fields[1]),
        qstart=int(fields[2]),
        qend=int(fields[3]),
        strand=strand,
        tname=fields[5],
        tlen=int(fields[6]),
        tstart=int(fields[7]),
        tend=int(fields[8]),
        nmatch=int(fields[9]),
        alen=int(fields[10]),
        mapq=int(fields[11]),
        is_primary=is_primary,
        nm=nm,
    )


def alignment_identity(record: PafRecord) -> float:
    """Per-alignment identity for gapit.reads/2: ``100 * (alen - nm) / alen``
    over the block length (column 10) and the NM tag (mismatches + gaps per
    the PAF spec). A row without NM cannot be assessed and counts as 100.0;
    a degenerate zero-length block also counts as 100.0 (no division)."""
    if record.nm is None or record.alen <= 0:
        return 100.0
    return 100.0 * (record.alen - record.nm) / record.alen


def filter_alignments(
    rows: Iterable[PafRecord], *, min_identity: float = 0.0, min_mapq: int = 0
) -> list[PafRecord]:
    """gapit.reads/2 row filter: keep alignments with identity >= min_identity
    and mapq >= min_mapq (applied after the primary-only rule upstream; both
    thresholds default off, in which case every row is kept verbatim — the
    gapit.reads/1 path must stay byte-identical, even for pathological rows
    whose arithmetic identity is negative)."""
    if min_identity <= 0.0 and min_mapq <= 0:
        return list(rows)
    return [row for row in rows if row.mapq >= min_mapq and alignment_identity(row) >= min_identity]


def union_length(intervals: Iterable[tuple[int, int]]) -> int:
    """Total length of the union of half-open [start, end) intervals:
    sort + linear sweep. Identical arithmetic to counting positions whose
    per-base depth is > 0, without materializing any per-base array."""
    covered = 0
    reach = -1
    for start, end in sorted(intervals):
        if end <= reach:
            continue
        covered += end - max(start, reach)
        reach = end
    return covered
