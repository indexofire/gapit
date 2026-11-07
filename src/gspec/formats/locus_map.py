"""Rich cluster locus-map JSON (prototype ``gspec.locus_map/1``).

This document is intentionally separate from gspec.cluster/1-2: compact
screen output stays byte-frozen, while callers that opt into detail get the
full LocusMap evidence.  A later cluster/3 can fold this into the screen
contract once the shape stabilizes.
"""

from __future__ import annotations

from collections.abc import Iterable
from datetime import datetime
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field

from gspec.cluster import BestCall, ClusterParams, ClusterReport
from gspec.formats.cluster import ClusterParamsDocument
from gspec.formats.json import ToolDocument, utc_timestamp
from gspec.locus_map_models import LocusMap
from gspec.schema_ids import Canonical
from gspec.typing_results import SchemeCall


class LocusMapFileDocument(BaseModel, frozen=True):
    """One file's detailed best-locus evidence."""

    file: str
    locus_maps: dict[str, LocusMap]
    phenotypes: dict[str, SchemeCall] = {}
    best: BestCall | None = None


class LocusMapDocument(BaseModel, frozen=True):
    """gspec.locus_map/1 — rich evidence for the best locus per file."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Annotated[Literal["gspec.locus_map/1"], Canonical] = Field(
        default="gspec.locus_map/1", alias="schema"
    )
    tool: ToolDocument = ToolDocument()
    created_at: str
    params: ClusterParamsDocument
    files: list[LocusMapFileDocument]


def render_locus_map_json(
    reports: Iterable[ClusterReport], params: ClusterParams, *, now: datetime
) -> str:
    """Serialize the opt-in rich evidence (indented, schema first)."""
    document = LocusMapDocument(
        created_at=utc_timestamp(now),
        params=ClusterParamsDocument(
            **params.model_dump(
                exclude={
                    "coverage_mode",
                    "fallback_evidence",
                    "locus_typeable_coverage",
                    "locus_typeable_identity",
                    "locus_quality_min_coverage",
                    "locus_quality_min_identity",
                    "locus_quality_max_missing_genes",
                    "locus_quality_max_extra_panel_genes",
                }
            )
        ),
        files=[
            LocusMapFileDocument(
                file=report.file,
                locus_maps=report.locus_maps,
                phenotypes=report.phenotypes,
                best=report.best,
            )
            for report in reports
            if report.locus_maps
        ],
    )
    return document.model_dump_json(indent=2, by_alias=True)
