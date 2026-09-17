"""FASTQ read screening via minimap2 (SPEC.md §10 — gapit extension)."""

import enum
import subprocess
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field

from gapit.db import Database
from gapit.dbcodec import decode_seqid
from gapit.errors import DependencyError, GapitError, InputError
from gapit.fasta import iter_fasta

ReadType = Literal["sr", "map-ont", "map-hifi"]


class ReadTypeEnum(enum.Enum):
    """CLI-facing read-type preset choices."""

    sr = "sr"
    map_ont = "map-ont"
    map_hifi = "map-hifi"


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


def _minimap2_version() -> str:
    """First line of ``minimap2 --version`` stdout ('' when minimap2 cannot
    run — the real invocation in run_minimap2 owns the typed error)."""
    try:
        result = subprocess.run(
            ["minimap2", "--version"], check=False, capture_output=True, text=True
        )
    except FileNotFoundError:
        return ""
    if result.returncode != 0:
        return ""
    return result.stdout.splitlines()[0].strip() if result.stdout else ""


def _mmi_index(database: Database) -> Path | None:
    """The persisted ``sequences.mmi`` to screen against, or None to use the
    FASTA.

    Usable iff the ``.mmi`` exists, a readable ``gapit-manifest.json`` sits
    beside it, and its ``minimap2_version`` matches the installed minimap2
    (index formats are version-specific — that equality is the compatibility
    gate; the sequences sha256 is deliberately NOT re-checked: the .mmi was
    built from the same ``sequences`` in the same directory). Any miss —
    missing file, unreadable/malformed manifest, empty or mismatched version
    — falls back to the FASTA silently: a performance fallback, not an error;
    genuine minimap2 failures surface in run_minimap2.
    """
    sequences = database.sequences_path
    mmi = sequences.with_name(f"{sequences.name}.mmi")
    if not mmi.is_file():
        return None
    # Imported here, not at module top: gapit.records -> formats.json ->
    # gapit.reads would be a circular import at load time.
    from gapit.records import read_manifest

    try:
        manifest = read_manifest(sequences.with_name("gapit-manifest.json"))
    except InputError:
        return None
    if not manifest.minimap2_version or manifest.minimap2_version != _minimap2_version():
        return None
    return mmi


def run_minimap2(
    lanes: list[tuple[Path, Path | None]],
    database: Database,
    *,
    read_type: ReadType,
    threads: int,
) -> list[PafRecord]:
    """Run one minimap2 invocation per lane (PAF on stdout) and concatenate
    the rows. The index argument is the persisted ``.mmi`` when usable (see
    ``_mmi_index``), else the FASTA (minimap2 then loads it per lane —
    negligible for the small dbs). minimap2's pairing semantics for >2 input
    files are undocumented; per-lane runs (r1[i] alone or with its mate
    r2[i]) are deterministic."""
    mmi = _mmi_index(database)
    index = database.sequences_path if mmi is None else mmi
    rows: list[PafRecord] = []
    for r1, r2 in lanes:
        argv = [
            "minimap2",
            "-x",
            read_type,
            "-t",
            str(threads),
            str(index),
            str(r1),
        ]
        if r2 is not None:
            argv.append(str(r2))
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
) -> ReadsReport:
    """Screen one sample's lanes against one database into a sample-level
    ReadsReport (union of all lanes' primary alignments)."""
    rows = run_minimap2(lanes, database, read_type=read_type, threads=threads)
    products = {record.id: record.description for record in iter_fasta(database.sequences_path)}
    genes = aggregate_coverage(
        rows, default_db=database.name, min_breadth=min_breadth, products=products
    )
    files = [str(r1) for r1, _ in lanes] + [str(r2) for _, r2 in lanes if r2 is not None]
    return ReadsReport(reads=tuple(files), genes=tuple(genes))
