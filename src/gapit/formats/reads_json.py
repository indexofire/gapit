"""Agent-facing reads JSON documents: gapit.reads/1, gapit.reads/2.

Split from formats/json.py at the 250-LOC ceiling: the reads-mode document
models and renderers (unfiltered reads/1 is the default and frozen; reads/2
adds the opt-in identity/MAPQ alignment filtering).
"""

from collections.abc import Iterable
from datetime import datetime
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from gapit.formats.json import ToolDocument, utc_timestamp
from gapit.reads import GeneCoverage, ReadsParams, ReadsReport


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
