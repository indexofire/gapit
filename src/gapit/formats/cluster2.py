"""gapit.cluster/2 renderers: multi-scheme cluster-screening output.

A cluster database whose typing document holds more than one scheme
(:func:`gapit.typing_models.uses_cluster_v2`) renders this VARIANT of the
cluster document: same body as /1, but each file's ``best`` drops its
phenotype fields and a ``phenotypes`` map keyed by scheme name (document
order) carries one :class:`~gapit.typing_results.SchemeCall` per scheme —
compose schemes included, rendered from their sibling calls. The TSV/CSV
header swaps PHENOTYPE for a packed PHENOTYPES column
(``o_group=O3;k_group=K6``; a phenotype containing ``;`` would degrade the
cell — the same accepted precedent as typing_result's GENES column).
Single-scheme and untyped databases keep byte-identical gapit.cluster/1
(formats/cluster.py). Registered in the schema registry as ``cluster2``.
"""

from collections.abc import Iterable
from datetime import datetime
from pathlib import PurePath
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from gapit.cluster import BestCall, ClusterParams, ClusterReport
from gapit.formats.cluster import (
    ClusterLocusDocument,
    ClusterParamsDocument,
    locus_document,
    missing_ids,
)
from gapit.formats.json import ToolDocument, utc_timestamp
from gapit.typing_results import SchemeCall

CLUSTER2_TSV_HEADER = (
    "FILE\tBEST_LOCUS\tTYPE\tPHENOTYPES\tCOVERAGE\tIDENTITY\tPRESENT\tPARTIAL\tMISSING_IDS"
)


class Cluster2BestDocument(BaseModel, frozen=True):
    """The best-locus summary without the phenotype fields — /2 moves every
    typing call into the file-level ``phenotypes`` map."""

    locus: str
    label: str
    type: str
    coverage_pct: float
    identity_pct: float
    genes_present: int
    genes_partial: int
    genes_absent: int


class Cluster2FileDocument(BaseModel, frozen=True):
    """One screened sample file: its best call (null below min_cluster_cov
    or with no locus coverage), the per-scheme typing calls in document
    order, and every covered locus, rank order."""

    file: str
    best: Cluster2BestDocument | None
    phenotypes: dict[str, SchemeCall] = Field(default_factory=dict)
    loci: list[ClusterLocusDocument]


class Cluster2Document(BaseModel, frozen=True):
    """gapit.cluster/2 — machine-readable multi-scheme cluster output."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.cluster/2"] = Field(default="gapit.cluster/2", alias="schema")
    tool: ToolDocument = ToolDocument()
    created_at: str
    params: ClusterParamsDocument
    files: list[Cluster2FileDocument]


def _best_document(best: BestCall) -> Cluster2BestDocument:
    return Cluster2BestDocument(
        locus=best.locus,
        label=best.label,
        type=best.type,
        coverage_pct=round(best.coverage_pct, 2),
        identity_pct=round(best.identity_pct, 2),
        genes_present=best.genes_present,
        genes_partial=best.genes_partial,
        genes_absent=best.genes_absent,
    )


def render_cluster2_json(
    reports: Iterable[ClusterReport], params: ClusterParams, *, now: datetime
) -> str:
    """Serialize cluster screening as gapit.cluster/2 (indented, schema
    first; ``best`` is null when no locus cleared min_cluster_cov — the
    phenotypes map still carries every scheme's call)."""
    document = Cluster2Document(
        created_at=utc_timestamp(now),
        params=ClusterParamsDocument(**params.model_dump()),
        files=[
            Cluster2FileDocument(
                file=report.file,
                best=None if report.best is None else _best_document(report.best),
                phenotypes=report.phenotypes,
                loci=[locus_document(locus) for locus in report.loci],
            )
            for report in reports
        ],
    )
    return document.model_dump_json(indent=2, by_alias=True)


def phenotypes_cell(report: ClusterReport) -> str:
    """The packed PHENOTYPES cell (the one rendering both the TSV column and
    the MD summary column consume): ``name=phenotype`` per scheme in
    document order, ``;``-joined; a null (ambiguous) call renders as
    ``name=-``."""
    return ";".join(f"{name}={call.phenotype or '-'}" for name, call in report.phenotypes.items())


def _tsv_row(report: ClusterReport, sep: str, nopath: bool) -> str:
    """One table line per file; a null best renders the /1 null-call dashes
    with the phenotypes cell still populated (schemes call without a best
    locus); missing ids join with ``;``."""
    file_column = PurePath(report.file).name if nopath else report.file
    phenotypes = phenotypes_cell(report)
    if report.best is None:
        return sep.join([file_column, "-", "-", phenotypes, "0.00", "0.00", "0", "0", "-"])
    best = report.best
    return sep.join(
        [
            file_column,
            best.locus,
            best.type,
            phenotypes,
            f"{best.coverage_pct:.2f}",
            f"{best.identity_pct:.2f}",
            str(best.genes_present),
            str(best.genes_partial),
            missing_ids(report, best.locus) or "-",
        ]
    )


def cluster2_tsv_preamble(*, csv: bool, noheader: bool) -> str:
    """The cluster/2 header chunk ("" under noheader; sinks skip empties)."""
    if noheader:
        return ""
    return CLUSTER2_TSV_HEADER.replace("\t", "," if csv else "\t") + "\n"


def cluster2_tsv_file_chunk(report: ClusterReport, *, csv: bool, nopath: bool) -> str:
    """One file's /2 chunk: its single table row (the streaming unit the
    cluster use-case emits as each file's report exists)."""
    return _tsv_row(report, "," if csv else "\t", nopath) + "\n"


def format_cluster2_tsv(
    reports: Iterable[ClusterReport],
    *,
    csv: bool,
    noheader: bool,
    nopath: bool,
) -> str:
    """Render cluster screening as TSV/CSV with the cluster/2 header (one
    row per file, PHENOTYPES packed per scheme)."""
    return cluster2_tsv_preamble(csv=csv, noheader=noheader) + "".join(
        cluster2_tsv_file_chunk(report, csv=csv, nopath=nopath) for report in reports
    )
