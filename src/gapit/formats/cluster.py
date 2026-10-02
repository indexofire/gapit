"""gapit.cluster/1 renderers: the cluster-screening JSON document and its
TSV/CSV table (stage 2 of the gene-cluster feature; the Markdown renderer
lives in formats/cluster_md.py, split at the 250-LOC ceiling in stage 3).

Typed databases (a ``typing.json`` beside ``sequences``) render a renderer
VARIANT: the TSV/CSV header gains a PHENOTYPE column after TYPE and the JSON
best block gains ``phenotype``/``phenotype_detail`` values — untyped output
stays byte-identical to the stage-2 documents (``phenotype_detail`` is
serialized only on typed runs). The schema registry (formats/schemas.py)
exposes the document as ``gapit schema cluster``.
"""

from collections.abc import Iterable
from datetime import datetime
from pathlib import PurePath
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from gapit.cluster import BestCall, ClusterParams, ClusterReport, GeneCall, LocusCall
from gapit.formats.json import ToolDocument, utc_timestamp
from gapit.typing_results import PhenotypeDetail

CLUSTER_TSV_HEADER = "FILE\tBEST_LOCUS\tTYPE\tCOVERAGE\tIDENTITY\tPRESENT\tPARTIAL\tMISSING_IDS"
CLUSTER_TSV_HEADER_TYPED = (
    "FILE\tBEST_LOCUS\tTYPE\tPHENOTYPE\tCOVERAGE\tIDENTITY\tPRESENT\tPARTIAL\tMISSING_IDS"
)


class ClusterGeneDocument(BaseModel, frozen=True):
    """One gene's verdict inside a screened locus."""

    gene_id: str
    start: int
    end: int
    strand: str
    coverage_pct: float
    identity_pct: float
    verdict: Literal["present", "partial", "absent"]


class ClusterLocusDocument(BaseModel, frozen=True):
    """One screened locus (rank order; only loci with coverage > 0)."""

    locus: str
    label: str
    type: str
    coverage_pct: float
    identity_pct: float
    rank: int
    genes: list[ClusterGeneDocument]
    missing: list[str]


class ClusterBestDocument(BaseModel, frozen=True):
    """The best-locus summary. ``phenotype``/``phenotype_detail`` carry the
    stage-3 typing call on typed databases (the detail key serializes only
    on typed runs; untyped documents keep the stage-2 byte shape)."""

    locus: str
    label: str
    type: str
    coverage_pct: float
    identity_pct: float
    genes_present: int
    genes_partial: int
    genes_absent: int
    phenotype: str | None = None
    phenotype_detail: PhenotypeDetail | None = None


class ClusterParamsDocument(BaseModel, frozen=True):
    """Cluster-screening parameters in effect (engine telemetry)."""

    db: str
    preset: str
    min_gene_cov: float
    min_gene_id: float
    min_cluster_cov: float
    threads: int


class ClusterFileDocument(BaseModel, frozen=True):
    """One screened sample file: its best call (null below min_cluster_cov
    or with no locus coverage) and every covered locus, rank order."""

    file: str
    best: ClusterBestDocument | None
    loci: list[ClusterLocusDocument]


class ClusterDocument(BaseModel, frozen=True):
    """gapit.cluster/1 — machine-readable cluster-screening output."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.cluster/1"] = Field(default="gapit.cluster/1", alias="schema")
    tool: ToolDocument = ToolDocument()
    created_at: str
    params: ClusterParamsDocument
    files: list[ClusterFileDocument]


def _gene_document(gene: GeneCall) -> ClusterGeneDocument:
    return ClusterGeneDocument(
        gene_id=gene.gene_id,
        start=gene.start,
        end=gene.end,
        strand=gene.strand,
        coverage_pct=round(gene.coverage_pct, 2),
        identity_pct=round(gene.identity_pct, 2),
        verdict=gene.verdict,
    )


def _locus_document(locus: LocusCall) -> ClusterLocusDocument:
    return ClusterLocusDocument(
        locus=locus.locus,
        label=locus.label,
        type=locus.type,
        coverage_pct=round(locus.coverage_pct, 2),
        identity_pct=round(locus.identity_pct, 2),
        rank=locus.rank,
        genes=[_gene_document(gene) for gene in locus.genes],
        missing=list(locus.missing),
    )


def _best_document(best: BestCall) -> ClusterBestDocument:
    return ClusterBestDocument(
        locus=best.locus,
        label=best.label,
        type=best.type,
        coverage_pct=round(best.coverage_pct, 2),
        identity_pct=round(best.identity_pct, 2),
        genes_present=best.genes_present,
        genes_partial=best.genes_partial,
        genes_absent=best.genes_absent,
        phenotype=best.phenotype,
        phenotype_detail=best.phenotype_detail,
    )


def render_cluster_json(
    reports: Iterable[ClusterReport], params: ClusterParams, *, now: datetime, typed: bool = False
) -> str:
    """Serialize cluster screening as gapit.cluster/1 (indented, schema
    first; ``best`` is null when no locus cleared min_cluster_cov). On
    untyped databases the additive ``phenotype_detail`` key is excluded so
    the document stays byte-identical to the stage-2 output; typed runs
    always carry it on every best block."""
    document = ClusterDocument(
        created_at=utc_timestamp(now),
        # field-for-field mirror of ClusterParams (engine telemetry contract)
        params=ClusterParamsDocument(**params.model_dump()),
        files=[
            ClusterFileDocument(
                file=report.file,
                best=None if report.best is None else _best_document(report.best),
                loci=[_locus_document(locus) for locus in report.loci],
            )
            for report in reports
        ],
    )
    if typed:
        return document.model_dump_json(indent=2, by_alias=True)
    return document.model_dump_json(
        indent=2,
        by_alias=True,
        exclude={"files": {"__all__": {"best": {"phenotype_detail"}}}},
    )


def _missing_ids(report: ClusterReport, locus_id: str) -> str:
    """The best locus's not-present gene ids, ``;``-joined (empty string
    when the call left nothing missing)."""
    return ";".join(
        gene_id for locus in report.loci if locus.locus == locus_id for gene_id in locus.missing
    )


def _tsv_row(report: ClusterReport, sep: str, nopath: bool, typed: bool) -> str:
    """One table line per file; a null best renders as the untyped row
    ``-  -  0.00  0.00  0  0  -`` (with a ``-`` phenotype on typed dbs);
    missing ids join with ``;``."""
    file_column = PurePath(report.file).name if nopath else report.file
    if report.best is None:
        cells = [file_column, "-", "-", "0.00", "0.00", "0", "0", "-"]
        if typed:
            cells.insert(3, "-")
        return sep.join(cells)
    best = report.best
    cells = [
        file_column,
        best.locus,
        best.type,
        f"{best.coverage_pct:.2f}",
        f"{best.identity_pct:.2f}",
        str(best.genes_present),
        str(best.genes_partial),
        _missing_ids(report, best.locus) or "-",
    ]
    if typed:
        cells.insert(3, best.phenotype or "-")
    return sep.join(cells)


def cluster_tsv_preamble(*, csv: bool, noheader: bool, typed: bool) -> str:
    """The cluster header chunk ("" under noheader; sinks skip empties)."""
    if noheader:
        return ""
    header = CLUSTER_TSV_HEADER_TYPED if typed else CLUSTER_TSV_HEADER
    return header.replace("\t", "," if csv else "\t") + "\n"


def cluster_tsv_file_chunk(report: ClusterReport, *, csv: bool, nopath: bool, typed: bool) -> str:
    """One file's chunk: its single table row (the streaming unit the
    cluster use-case emits as each file's report exists)."""
    return _tsv_row(report, "," if csv else "\t", nopath, typed) + "\n"


def format_cluster_tsv(
    reports: Iterable[ClusterReport],
    *,
    csv: bool,
    noheader: bool,
    nopath: bool,
    typed: bool = False,
) -> str:
    """Render cluster screening as TSV/CSV with the cluster-engine header
    (one row per file; abricate-unrelated by design). Typed databases gain
    the PHENOTYPE column after TYPE; untyped output keeps the stage-2
    header and rows byte-identically."""
    return cluster_tsv_preamble(csv=csv, noheader=noheader, typed=typed) + "".join(
        cluster_tsv_file_chunk(report, csv=csv, nopath=nopath, typed=typed) for report in reports
    )
