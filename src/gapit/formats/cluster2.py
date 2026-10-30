"""gapit.cluster/2 renderers: multi-scheme cluster-screening output.

A cluster database whose typing document holds more than one scheme
(:func:`gapit.typing_models.uses_cluster_v2`) renders this VARIANT of the
cluster document: same body as /1, but each file's ``best`` drops its
phenotype fields and a ``phenotypes`` map keyed by scheme name (document
order) carries one :class:`~gapit.typing_results.SchemeCall` per scheme —
compose schemes included, rendered from their sibling calls, each
rule-based call carrying its winning locus's evidence (the additive
``winner_*`` fields).

The TSV/CSV table is DERIVED from the typing document (no global-best
columns): ``FILE``, then every compose scheme's value column in document
order (uppercased name, ``-``→``_``), then each non-compose scheme's
seven-column block in document order — the call, then the winning
locus's evidence (``<P>_LOCUS``, ``<P>_COVERAGE``, ``<P>_IDENTITY``,
``<P>_PRESENT``, ``<P>_PARTIAL``, ``<P>_MISSING`` where ``<P>`` is the
uppercased scheme name minus a trailing ``_GROUP``, so ``o_group``
renders ``O_GROUP`` then ``O_LOCUS``). Fallback calls render the
fallback string with ``-`` evidence; ambiguous calls render ``-`` in the
call column. For vpa_ok the header is ``FILE  SEROTYPE  O_GROUP  O_LOCUS
… K_MISSING`` — the composed serotype leads. Single-scheme and untyped
databases keep byte-identical gapit.cluster/1 (formats/cluster.py).
Registered in the schema registry as ``cluster2``.
"""

from collections import Counter
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
)
from gapit.formats.csvutil import delimited_line
from gapit.formats.json import ToolDocument, utc_timestamp
from gapit.typing_models import TypingDocument
from gapit.typing_results import SchemeCall

_EVIDENCE_COLUMNS = ("LOCUS", "COVERAGE", "IDENTITY", "PRESENT", "PARTIAL", "MISSING")
_GROUP_SUFFIX = "_GROUP"


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


def _scheme_column(name: str) -> str:
    """A scheme name as a column token: uppercased, ``-`` → ``_``."""
    return name.upper().replace("-", "_")


def _short_evidence_prefix(name: str) -> str:
    """The preferred evidence prefix: the token minus a trailing ``_GROUP``."""
    token = _scheme_column(name)
    return (
        token[: -len(_GROUP_SUFFIX)]
        if token.endswith(_GROUP_SUFFIX) and len(token) > len(_GROUP_SUFFIX)
        else token
    )


def _evidence_prefix(name: str, *, colliding: frozenset[str] = frozenset()) -> str:
    """The evidence-column prefix, disambiguated when short tokens collide.

    Two schemes such as ``o_group`` and ``o`` both prefer ``O``; in that case
    both keep their full token (``O_GROUP`` / ``O``) so every column stays
    unique. The common non-colliding case keeps the familiar short form.
    """
    token = _scheme_column(name)
    short = _short_evidence_prefix(name)
    return token if short in colliding else short


def cluster2_columns(document: TypingDocument) -> tuple[str, ...]:
    """The /2 TSV/CSV column names: ``FILE``, every compose scheme's value
    column in document order, then each non-compose scheme's call +
    winner-evidence block in document order (deterministic in the doc)."""
    columns: list[str] = ["FILE"]
    columns.extend(
        _scheme_column(scheme.name) for scheme in document.schemes if scheme.compose is not None
    )
    non_compose = [scheme.name for scheme in document.schemes if scheme.compose is None]
    short_counts = Counter(_short_evidence_prefix(name) for name in non_compose)
    colliding = frozenset(short for short, count in short_counts.items() if count > 1)
    for scheme in document.schemes:
        if scheme.compose is None:
            prefix = _evidence_prefix(scheme.name, colliding=colliding)
            columns.append(_scheme_column(scheme.name))
            columns.extend(f"{prefix}_{suffix}" for suffix in _EVIDENCE_COLUMNS)
    return tuple(columns)


def pct_or_dash(value: float | None) -> str:
    """A percentage cell: ``%.2f`` or ``-`` (the shared winner-evidence
    formatting of the TSV and MD /2 renderers)."""
    return "-" if value is None else f"{value:.2f}"


def count_or_dash(value: int | None) -> str:
    """A count cell: the integer or ``-`` (unset evidence)."""
    return "-" if value is None else str(value)


def cluster2_row_cells(
    report: ClusterReport, document: TypingDocument, *, nopath: bool
) -> list[str]:
    """One file's /2 cells aligned to :func:`cluster2_columns` (the shared
    builder behind the TSV row and the MD summary row): compose values,
    then each scheme's call followed by its winning locus's evidence —
    ``-`` cells when the call is fallback/gate-based or ambiguous."""
    cells = [PurePath(report.file).name if nopath else report.file]
    cells.extend(
        report.phenotypes[scheme.name].phenotype or "-"
        for scheme in document.schemes
        if scheme.compose is not None
    )
    for scheme in document.schemes:
        if scheme.compose is None:
            call = report.phenotypes[scheme.name]
            cells.extend(
                (
                    call.phenotype or "-",
                    call.winner_locus or "-",
                    pct_or_dash(call.winner_coverage_pct),
                    pct_or_dash(call.winner_identity_pct),
                    count_or_dash(call.winner_present),
                    count_or_dash(call.winner_partial),
                    ";".join(call.winner_missing) if call.winner_missing else "-",
                )
            )
    return cells


def cluster2_tsv_preamble(*, csv: bool, noheader: bool, document: TypingDocument) -> str:
    """The cluster/2 header chunk ("" under noheader; sinks skip empties)."""
    if noheader:
        return ""
    return delimited_line(cluster2_columns(document), as_csv=csv)


def cluster2_tsv_file_chunk(
    report: ClusterReport, *, csv: bool, nopath: bool, document: TypingDocument
) -> str:
    """One file's /2 chunk: its single table row (the streaming unit the
    cluster use-case emits as each file's report exists)."""
    return delimited_line(cluster2_row_cells(report, document, nopath=nopath), as_csv=csv)


def format_cluster2_tsv(
    reports: Iterable[ClusterReport],
    *,
    csv: bool,
    noheader: bool,
    nopath: bool,
    document: TypingDocument,
) -> str:
    """Render cluster screening as TSV/CSV with the cluster/2 header (one
    row per file: compose schemes first, then each scheme's call and its
    winning-locus evidence)."""
    return cluster2_tsv_preamble(csv=csv, noheader=noheader, document=document) + "".join(
        cluster2_tsv_file_chunk(report, csv=csv, nopath=nopath, document=document)
        for report in reports
    )
