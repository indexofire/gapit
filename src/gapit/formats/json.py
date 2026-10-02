"""Agent-facing JSON documents: gapit.report/1, gapit.version/1.

The reads-mode documents (gapit.reads/1 and /2) live in
formats/reads_json.py — split at the 250-LOC ceiling — and the typing
designation document (gapit.typing_result/1) in formats/typing_result.py.
gapit.report/1 is pure gene detection: typed and untyped gene databases
serialize byte-identically.
"""

from collections.abc import Iterable
from datetime import UTC, datetime
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from gapit import __version__
from gapit.hits import Hit
from gapit.report import Report, ScreeningParams


class FragmentDocument(BaseModel, frozen=True):
    """One contributing fragment of a merged hit (merge mode only)."""

    contig: str
    start: int
    end: int
    strand: str
    identity_pct: float
    coverage_pct: float


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
    # Additive gapit extension (schema-compatible, no version bump): present
    # only on --merge-fragments rows; None is excluded at serialization, so
    # default-mode documents stay byte-identical.
    merged: bool | None = None
    fragments: list[FragmentDocument] | None = None


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


class VersionDocument(BaseModel, frozen=True):
    """gapit.version/1 — compact one-line self-description."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.version/1"] = Field(default="gapit.version/1", alias="schema")
    name: str = "gapit"
    version: str


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
        merged=True if hit.merged else None,
        fragments=(
            [
                FragmentDocument(
                    contig=fragment.contig,
                    start=fragment.start,
                    end=fragment.end,
                    strand=fragment.strand,
                    identity_pct=round(fragment.identity_pct, 2),
                    coverage_pct=round(fragment.coverage_pct, 2),
                )
                for fragment in hit.fragments
            ]
            or None
        ),
    )


def render_json(reports: Iterable[Report], params: ScreeningParams, *, now: datetime) -> str:
    """Serialize screening results as gapit.report/1 (indented, schema first)."""
    document = ReportDocument(
        created_at=utc_timestamp(now),
        params=ParamsDocument(
            db=params.db, minid=params.minid, mincov=params.mincov, threads=params.threads
        ),
        files=[
            FileDocument(
                file=report.file,
                hits=[_hit_document(hit) for hit in report.hits],
            )
            for report in reports
        ],
    )
    return document.model_dump_json(indent=2, by_alias=True, exclude_none=True)
