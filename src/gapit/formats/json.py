"""Agent-facing JSON documents: gapit.report/1, gapit.list/1, gapit.version/1."""

from collections.abc import Iterable
from datetime import UTC, datetime
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from gapit import __version__
from gapit.hits import Hit
from gapit.reads import GeneCoverage, ReadsParams, ReadsReport
from gapit.report import Report, ScreeningParams


class HitDocument(BaseModel, frozen=True):
    """One hit — 1:1 with the TSV columns, snake_case, same string values."""

    sequence: str
    start: int
    end: int
    strand: str
    gene: str
    coverage: str
    coverage_map: str
    gaps: str
    coverage_pct: float
    identity_pct: float
    database: str
    accession: str
    product: str
    # "resistance" is frozen by gapit.report/1 (1:1 with the TSV RESISTANCE
    # column); the value flows from Hit.function and carries functional
    # categories for native DBs (Wave F1 renamed the internal slot only).
    resistance: str


class FileDocument(BaseModel, frozen=True):
    """One screened file and its hits (zero hits -> empty list)."""

    file: str
    hits: list[HitDocument]


class ParamsDocument(BaseModel, frozen=True):
    """The screening parameters in effect."""

    db: str
    minid: float
    mincov: float
    threads: int


class ToolDocument(BaseModel, frozen=True):
    """Tool self-identification."""

    name: str = "gapit"
    version: str = __version__


class ReportDocument(BaseModel, frozen=True):
    """gapit.report/1 — the canonical machine-readable screening output."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.report/1"] = Field(default="gapit.report/1", alias="schema")
    tool: ToolDocument = ToolDocument()
    created_at: str
    params: ParamsDocument
    files: list[FileDocument]


class ListEntryDocument(BaseModel, frozen=True):
    """One row of `gapit list --json`."""

    name: str
    sequences: int
    dbtype: str
    date: str


class ListDocument(BaseModel, frozen=True):
    """gapit.list/1."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.list/1"] = Field(default="gapit.list/1", alias="schema")
    databases: list[ListEntryDocument]


class VersionDocument(BaseModel, frozen=True):
    """gapit.version/1 — compact one-line self-description."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.version/1"] = Field(default="gapit.version/1", alias="schema")
    name: str = "gapit"
    version: str


class GeneCoverageDocument(BaseModel, frozen=True):
    """One gene's presence call in reads mode."""

    gene: str
    database: str
    accession: str
    product: str
    # "resistance" is frozen by gapit.reads/1; the value flows from
    # GeneCoverage.function (functional categories for native DBs).
    resistance: str
    tlen: int
    breadth_pct: float
    mean_depth: float
    reads_mapped: int
    present: bool


class ReadsFileDocument(BaseModel, frozen=True):
    """One screened read set."""

    reads: list[str]
    genes: list[GeneCoverageDocument]


class ReadsParamsDocument(BaseModel, frozen=True):
    """Read-screening parameters in effect."""

    db: str
    read_type: str
    min_breadth: float
    threads: int


class ReadsDocument(BaseModel, frozen=True):
    """gapit.reads/1 — machine-readable read-screening output."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.reads/1"] = Field(default="gapit.reads/1", alias="schema")
    tool: ToolDocument = ToolDocument()
    created_at: str
    params: ReadsParamsDocument
    files: list[ReadsFileDocument]


class GeneCoverage2Document(BaseModel, frozen=True):
    """One gene's presence call in reads/2 mode: the reads/1 fields plus the
    alen-weighted mean per-alignment identity."""

    gene: str
    database: str
    accession: str
    product: str
    resistance: str
    tlen: int
    breadth_pct: float
    mean_depth: float
    reads_mapped: int
    present: bool
    mean_identity_pct: float


class Reads2FileDocument(BaseModel, frozen=True):
    """One screened read set (reads/2)."""

    reads: list[str]
    genes: list[GeneCoverage2Document]


class Reads2ParamsDocument(BaseModel, frozen=True):
    """Read-screening parameters in effect (reads/2): reads/1 params plus the
    opt-in alignment filters."""

    db: str
    read_type: str
    min_breadth: float
    threads: int
    min_identity: float
    min_mapq: int


class Reads2Document(BaseModel, frozen=True):
    """gapit.reads/2 — read-screening output under opt-in identity/MAPQ
    alignment filtering (reads/1 stays the default and is frozen)."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.reads/2"] = Field(default="gapit.reads/2", alias="schema")
    tool: ToolDocument = ToolDocument()
    created_at: str
    params: Reads2ParamsDocument
    files: list[Reads2FileDocument]


def utc_timestamp(now: datetime) -> str:
    """ISO-8601 UTC with a trailing Z, second precision."""
    return now.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _hit_document(hit: Hit) -> HitDocument:
    return HitDocument(
        sequence=hit.sequence,
        start=hit.start,
        end=hit.end,
        strand=hit.strand,
        gene=hit.gene,
        coverage=f"{hit.s_start}-{hit.s_end}/{hit.s_len}",
        coverage_map=hit.coverage_map,
        gaps=f"{hit.gap_openings}/{hit.gaps}",
        coverage_pct=round(hit.coverage_pct, 2),
        identity_pct=round(hit.identity_pct, 2),
        database=hit.database,
        accession=hit.accession,
        product=hit.product,
        resistance=hit.function,
    )


def render_json(reports: Iterable[Report], params: ScreeningParams, *, now: datetime) -> str:
    """Serialize screening results as gapit.report/1 (indented, schema first)."""
    document = ReportDocument(
        created_at=utc_timestamp(now),
        params=ParamsDocument(
            db=params.db, minid=params.minid, mincov=params.mincov, threads=params.threads
        ),
        files=[
            FileDocument(file=report.file, hits=[_hit_document(hit) for hit in report.hits])
            for report in reports
        ],
    )
    return document.model_dump_json(indent=2, by_alias=True)


def _gene_coverage_document(gene: GeneCoverage) -> GeneCoverageDocument:
    return GeneCoverageDocument(
        gene=gene.gene,
        database=gene.database,
        accession=gene.accession,
        product=gene.product,
        resistance=gene.function,
        tlen=gene.tlen,
        breadth_pct=round(gene.breadth_pct, 2),
        mean_depth=round(gene.mean_depth, 2),
        reads_mapped=gene.reads_mapped,
        present=gene.present,
    )


def render_reads_json(reports: Iterable[ReadsReport], params: ReadsParams, *, now: datetime) -> str:
    """Serialize read-screening results as gapit.reads/1 (indented, schema first)."""
    document = ReadsDocument(
        created_at=utc_timestamp(now),
        params=ReadsParamsDocument(
            db=params.db,
            read_type=params.read_type,
            min_breadth=params.min_breadth,
            threads=params.threads,
        ),
        files=[
            ReadsFileDocument(
                reads=list(report.reads),
                genes=[_gene_coverage_document(gene) for gene in report.genes],
            )
            for report in reports
        ],
    )
    return document.model_dump_json(indent=2, by_alias=True)


def _gene_coverage2_document(gene: GeneCoverage) -> GeneCoverage2Document:
    return GeneCoverage2Document(
        gene=gene.gene,
        database=gene.database,
        accession=gene.accession,
        product=gene.product,
        resistance=gene.function,
        tlen=gene.tlen,
        breadth_pct=round(gene.breadth_pct, 2),
        mean_depth=round(gene.mean_depth, 2),
        reads_mapped=gene.reads_mapped,
        present=gene.present,
        mean_identity_pct=round(gene.mean_identity_pct, 2),
    )


def render_reads2_json(
    reports: Iterable[ReadsReport], params: ReadsParams, *, now: datetime
) -> str:
    """Serialize filtered read-screening results as gapit.reads/2 (indented,
    schema first; same shape as /1 plus the filter params and per-gene
    mean_identity_pct)."""
    document = Reads2Document(
        created_at=utc_timestamp(now),
        params=Reads2ParamsDocument(
            db=params.db,
            read_type=params.read_type,
            min_breadth=params.min_breadth,
            threads=params.threads,
            min_identity=params.min_identity,
            min_mapq=params.min_mapq,
        ),
        files=[
            Reads2FileDocument(
                reads=list(report.reads),
                genes=[_gene_coverage2_document(gene) for gene in report.genes],
            )
            for report in reports
        ],
    )
    return document.model_dump_json(indent=2, by_alias=True)
