"""FASTQ/FASTA read screening via minimap2 (SPEC.md §10 — gapit extension)."""

import enum
import gzip
import shlex
import subprocess
import sys
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Literal, Protocol

from pydantic import BaseModel, Field

from gapit.db import Database
from gapit.dbcodec import decode_seqid
from gapit.errors import DependencyError, GapitError, InputError
from gapit.fasta import iter_fasta

ReadType = Literal["sr", "map-ont", "map-hifi"]

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
    natively, so detection must not reject them."""
    with path.open("rb") as raw:
        compressed = raw.read(2) == _GZIP_MAGIC
    if compressed:
        with gzip.open(path, "rb") as handle:
            return _peek_read_kind(handle, path)
    with path.open("rb") as handle:
        return _peek_read_kind(handle, path)


class ReadsParams(BaseModel, frozen=True):
    """Parameters for one read-screening run (SPEC.md §10)."""

    db: str
    read_type: ReadType = "sr"
    min_breadth: float = Field(default=90.0, ge=0.0, le=100.0)
    threads: int = Field(default=1, ge=1)


class PafRecord(BaseModel, frozen=True):
    """One PAF alignment row: 12 required fields + primary flag (tp:A:P)."""

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


def parse_paf_row(line: str) -> PafRecord:
    """Parse one tab-delimited PAF line; <12 fields is a hard error. Records
    without a tp tag count as primary."""
    fields = line.split("\t")
    if len(fields) < 12:
        raise GapitError(
            f"can not parse PAF row (expected >= 12 fields): {line!r}",
            code="PAF_PARSE_FAILED",
        )
    is_primary = True
    for tag in fields[12:]:
        if tag.startswith("tp:A:"):
            is_primary = tag == "tp:A:P"
    strand = fields[4]
    if strand not in ("+", "-"):
        raise GapitError(
            f"can not parse PAF strand (expected + or -): {strand!r}",
            code="PAF_PARSE_FAILED",
        )
    return PafRecord(
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
    )


class GeneCoverage(BaseModel, frozen=True):
    """Per-gene presence call over primary alignments."""

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
) -> list[GeneCoverage]:
    """Aggregate primary alignments per target into exact per-base coverage:
    breadth = union of tstart..tend, depth = per-base sum, reads = distinct
    qnames. Zero-read genes are omitted; output sorts by (-breadth, gene)."""
    depths: dict[str, list[int]] = {}
    qnames: dict[str, set[str]] = {}
    for row in rows:
        if not row.is_primary:
            continue
        if row.tname not in depths:
            depths[row.tname] = [0] * row.tlen
            qnames[row.tname] = set()
        depth = depths[row.tname]
        for position in range(row.tstart, row.tend):
            depth[position] += 1
        qnames[row.tname].add(row.qname)
    descriptions = products or {}
    coverages: list[GeneCoverage] = []
    for tname, depth in depths.items():
        tlen = len(depth)
        covered = sum(1 for value in depth if value > 0)
        header = decode_seqid(tname, default_db)
        coverages.append(
            GeneCoverage(
                database=header.database,
                gene=header.gene,
                accession=header.accession,
                function=header.function,
                product=descriptions.get(tname, ""),
                tlen=tlen,
                breadth_pct=100.0 * covered / tlen,
                mean_depth=sum(depth) / tlen,
                reads_mapped=len(qnames[tname]),
                present=100.0 * covered / tlen >= min_breadth,
            )
        )
    return sorted(coverages, key=lambda entry: (-entry.breadth_pct, entry.gene))


def run_minimap2(
    lanes: list[tuple[Path, Path | None]],
    database: Database,
    *,
    read_type: ReadType,
    threads: int,
    debug: bool = False,
) -> list[PafRecord]:
    """Run one minimap2 invocation per lane (PAF on stdout) and concatenate
    the rows. The index argument is always the ``sequences`` FASTA — minimap2
    indexes it in memory with the invocation preset's own parameters. A
    persisted ``.mmi`` is deliberately rejected even when one sits beside the
    FASTA: a default-built index overrides the ``-x`` preset's indexing
    parameters (``-k, -w or -H overridden by prebuilt index``), which
    misassigns close homologs and benchmarks slower than in-memory indexing
    (2026-09-19: blaCTX-M/blaSHV allele divergence, +1.2 s on the ncbi db).
    minimap2's pairing semantics for >2 input files are undocumented;
    per-lane runs (r1[i] alone or with its mate r2[i]) are deterministic.
    With ``debug``, echo each argv to stderr (abricate --debug parity)."""
    rows: list[PafRecord] = []
    for r1, r2 in lanes:
        argv = [
            "minimap2",
            "-x",
            read_type,
            "-t",
            str(threads),
            str(database.sequences_path),
            str(r1),
        ]
        if r2 is not None:
            argv.append(str(r2))
        if debug:
            print(f"gapit: run: {shlex.join(argv)}", file=sys.stderr)
        try:
            result = subprocess.run(argv, check=False, capture_output=True, text=True)
        except FileNotFoundError as exc:
            raise DependencyError(
                "required binary not found on PATH: minimap2",
                code="MISSING_DEPENDENCY",
                context={"binary": "minimap2"},
            ) from exc
        if result.returncode != 0:
            raise GapitError(
                f"minimap2 failed: {result.stderr.strip()}",
                code="MINIMAP2_FAILED",
                context={"binary": "minimap2", "file": str(r1)},
            )
        rows.extend(parse_paf_row(line) for line in result.stdout.splitlines() if line.strip())
    return rows


def screen_reads(
    lanes: list[tuple[Path, Path | None]],
    database: Database,
    *,
    read_type: ReadType,
    min_breadth: float,
    threads: int,
    debug: bool = False,
) -> ReadsReport:
    """Screen one sample's lanes against one database into a sample-level
    ReadsReport (union of all lanes' primary alignments)."""
    rows = run_minimap2(lanes, database, read_type=read_type, threads=threads, debug=debug)
    products = {record.id: record.description for record in iter_fasta(database.sequences_path)}
    genes = aggregate_coverage(
        rows, default_db=database.name, min_breadth=min_breadth, products=products
    )
    files = [str(r1) for r1, _ in lanes] + [str(r2) for _, r2 in lanes if r2 is not None]
    return ReadsReport(reads=tuple(files), genes=tuple(genes))
