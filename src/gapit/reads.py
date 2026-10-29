"""FASTQ/FASTA read screening via minimap2 (SPEC.md §10 — gapit extension)."""

import enum
import gzip
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Protocol

from pydantic import BaseModel, Field

from gapit.db import Database
from gapit.dbcodec import decode_seqid
from gapit.errors import InputError
from gapit.fasta import iter_fasta_headers
from gapit.gene_floors import gene_floors, iter_apply_gene_floors
from gapit.minimap2_run import iter_minimap2
from gapit.paf import (
    PafRecord,
    ReadType,
    alignment_identity,
    iter_filter_alignments,
    union_length,
)

_GZIP_MAGIC = b"\x1f\x8b"


class _ByteStream(Protocol):
    """The byte-read surface shared by raw and gzip-decompressed handles."""

    def read(self, size: int = -1, /) -> bytes: ...


class ReadTypeEnum(enum.Enum):
    """CLI-facing read-type preset choices."""

    sr = "sr"
    map_ont = "map-ont"
    map_hifi = "map-hifi"


class ReadFileKind(enum.Enum):
    """Input file kind detected from content (minimap2 takes FASTA and FASTQ
    queries natively; gapit detects to resolve the preset)."""

    fasta = "fasta"
    fastq = "fastq"


def _peek_read_kind(handle: _ByteStream, path: Path) -> ReadFileKind:
    """First non-whitespace byte: '>' = FASTA, '@' = FASTQ; anything else (or
    EOF) is a typed input error."""
    while True:
        byte = handle.read(1)
        if not byte:
            raise InputError(
                f"reads file is empty: {path}",
                code="INVALID_READS_FORMAT",
                context={"file": str(path)},
            )
        if byte.isspace():
            continue
        if byte == b">":
            return ReadFileKind.fasta
        if byte == b"@":
            return ReadFileKind.fastq
        raise InputError(
            f"reads file is neither FASTA nor FASTQ: {path}",
            code="INVALID_READS_FORMAT",
            context={"file": str(path)},
        )


def detect_read_kind(path: Path) -> ReadFileKind:
    """Detect a --r1/--r2 file's kind from content. Gzip-wrapped files
    (magic 1f 8b) are peeked through the decompressor: minimap2 reads them
    natively, so detection must not reject them. A truncated or corrupt
    gzip stream raises the same typed INVALID_READS_FORMAT error as any
    other unreadable reads file, never a raw EOFError."""
    with path.open("rb") as raw:
        compressed = raw.read(2) == _GZIP_MAGIC
    if compressed:
        with gzip.open(path, "rb") as handle:
            try:
                return _peek_read_kind(handle, path)
            except (OSError, EOFError, UnicodeDecodeError) as exc:
                raise InputError(
                    f"cannot decompress reads file: {path}",
                    code="INVALID_READS_FORMAT",
                    context={"file": str(path)},
                ) from exc
    with path.open("rb") as handle:
        return _peek_read_kind(handle, path)


class ReadsParams(BaseModel, frozen=True):
    """Parameters for one read-screening run (SPEC.md §10). The reads/2
    thresholds default off; gapit.reads/1 rendering ignores them entirely.
    ``min_depth``/``min_relative_depth`` add the abundance gate: absolute
    mean depth and the gene:control mean-depth ratio (0..1)."""

    db: str
    read_type: ReadType = "sr"
    min_breadth: float = Field(default=90.0, ge=0.0, le=100.0)
    threads: int = Field(default=1, ge=1)
    min_identity: float = Field(default=0.0, ge=0.0, le=100.0)
    min_mapq: int = Field(default=0, ge=0)
    min_depth: float = Field(default=0.0, ge=0.0)
    min_relative_depth: float = Field(default=0.0, ge=0.0, le=1.0)
    relative_to: str = "lpxA"


class GeneCoverage(BaseModel, frozen=True):
    """Per-gene presence call over primary alignments. mean_identity_pct is
    the alen-weighted mean of per-alignment identity over the aggregated
    (kept) rows — the weighting favors long alignments; computed in reads/1
    mode too but rendered only by gapit.reads/2. ``relative_depth_pct`` is
    the gene:control mean-depth ratio as a percentage (None when no control
    gene was found or relative-depth gating is off)."""

    database: str
    gene: str
    accession: str
    function: str
    product: str
    tlen: int
    breadth_pct: float
    mean_depth: float
    reads_mapped: int
    present: bool
    mean_identity_pct: float
    relative_depth_pct: float | None = None


class ReadsReport(BaseModel, frozen=True):
    """Result of screening one set of read files."""

    reads: tuple[str, ...]
    genes: tuple[GeneCoverage, ...]


def aggregate_coverage(
    rows: Iterable[PafRecord],
    *,
    default_db: str,
    min_breadth: float,
    products: Mapping[str, str] | None = None,
    min_depth: float = 0.0,
    min_relative_depth: float = 0.0,
    relative_to: str | None = None,
) -> list[GeneCoverage]:
    """Aggregate primary alignments per target: breadth = union of
    tstart..tend intervals, depth = summed interval lengths, reads = distinct
    qnames. tlen is the FIRST-SEEN row's tlen per target (denominator
    contract). Zero-read genes are omitted; output sorts by (-breadth, gene).

    ``min_depth``/``min_relative_depth`` gate presence after aggregation.
    The relative gate compares each gene's mean depth with the mean depth of
    ``relative_to`` (default lpxA at the call sites); a missing control with
    ``min_relative_depth > 0`` makes every call absent (the control is the
    sample-level gate, so its absence must not silently pass everything).
    """
    intervals: dict[str, list[tuple[int, int]]] = {}
    tlens: dict[str, int] = {}
    qnames: dict[str, set[str]] = {}
    weights: dict[str, list[tuple[float, int]]] = {}
    for row in rows:
        if not row.is_primary:
            continue
        if row.tname not in intervals:
            intervals[row.tname] = []
            tlens[row.tname] = row.tlen
            qnames[row.tname] = set()
            weights[row.tname] = []
        intervals[row.tname].append((row.tstart, row.tend))
        qnames[row.tname].add(row.qname)
        weights[row.tname].append((alignment_identity(row), row.alen))
    descriptions = products or {}
    coverages: list[GeneCoverage] = []
    for tname, spans in intervals.items():
        tlen = tlens[tname]
        covered = union_length(spans)
        breadth_pct = 100.0 * covered / tlen
        header = decode_seqid(tname, default_db)
        pairs = weights[tname]
        weight_sum = sum(weight for _, weight in pairs)
        coverages.append(
            GeneCoverage(
                database=header.database,
                gene=header.gene,
                accession=header.accession,
                function=header.function,
                product=descriptions.get(tname, ""),
                tlen=tlen,
                breadth_pct=breadth_pct,
                mean_depth=sum(end - start for start, end in spans) / tlen,
                reads_mapped=len(qnames[tname]),
                present=breadth_pct >= min_breadth,
                mean_identity_pct=(
                    sum(identity * weight for identity, weight in pairs) / weight_sum
                    if weight_sum
                    else 0.0
                ),
            )
        )
    control_depth = (
        max((entry.mean_depth for entry in coverages if entry.gene == relative_to), default=None)
        if relative_to
        else None
    )
    gated: list[GeneCoverage] = []
    for entry in coverages:
        relative_pct = (
            None if control_depth in (None, 0.0) else 100.0 * entry.mean_depth / control_depth
        )
        present = entry.present and entry.mean_depth >= min_depth
        if min_relative_depth > 0.0:
            present = (
                present and relative_pct is not None and relative_pct >= 100.0 * min_relative_depth
            )
        gated.append(
            entry.model_copy(update={"present": present, "relative_depth_pct": relative_pct})
        )
    return sorted(gated, key=lambda entry: (-entry.breadth_pct, entry.gene))


def screen_reads(
    lanes: list[tuple[Path, Path | None]],
    database: Database,
    *,
    read_type: ReadType,
    min_breadth: float,
    threads: int,
    debug: bool = False,
    min_identity: float = 0.0,
    min_mapq: int = 0,
    min_depth: float = 0.0,
    min_relative_depth: float = 0.0,
    relative_to: str = "lpxA",
) -> ReadsReport:
    """Screen one sample's lanes against one database into a sample-level
    ReadsReport (union of all lanes' primary alignments). With a nonzero
    min_identity/min_mapq (gapit.reads/2), alignments are filtered BEFORE
    aggregation and the identity rule needs the NM tags. ``nm_tags`` is
    ALWAYS on: minimap2's no-CIGAR mode emits clipped coordinates for
    diverged alignments (26ECO0071 astA: 57% breadth without --cs, 100%
    with), silently under-reporting breadth — the reads/1 core metric.

    A database carrying a gapit.floors/1 sidecar additionally drops
    sub-floor alignments per gene BEFORE aggregation
    (:func:`gapit.gene_floors.apply_gene_floors`) — db-driven opt-in; a
    floorless database is byte-identical to the pre-floors behavior."""
    reads2 = min_identity > 0.0 or min_mapq > 0
    rows: Iterable[PafRecord] = iter_minimap2(
        lanes,
        database,
        read_type=read_type,
        threads=threads,
        debug=debug,
        nm_tags=True,
    )
    if reads2:
        rows = iter_filter_alignments(rows, min_identity=min_identity, min_mapq=min_mapq)
    floors = gene_floors(database)
    if floors is not None:
        rows = iter_apply_gene_floors(rows, floors, default_db=database.name)
    products = dict(iter_fasta_headers(database.sequences_path))
    genes = aggregate_coverage(
        rows,
        default_db=database.name,
        min_breadth=min_breadth,
        products=products,
        min_depth=min_depth,
        min_relative_depth=min_relative_depth,
        relative_to=relative_to,
    )
    files = [str(r1) for r1, _ in lanes] + [str(r2) for _, r2 in lanes if r2 is not None]
    return ReadsReport(reads=tuple(files), genes=tuple(genes))
