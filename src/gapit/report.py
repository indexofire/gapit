"""Screening parameters and the canonical Report model."""

from pydantic import BaseModel, Field

from gapit.hits import Hit


class ScreeningParams(BaseModel, frozen=True):
    """Parameters for one screening run (SPEC.md §1 defaults and bounds)."""

    db: str
    minid: float = Field(default=80.0, gt=0.0, le=100.0)
    mincov: float = Field(default=80.0, ge=0.0, le=100.0)
    threads: int = Field(default=1, ge=1)


class Report(BaseModel, frozen=True):
    """Canonical in-memory result of screening one input file.

    Hits are sorted (sequence, start) before construction; the tuple is the
    final, stable order.
    """

    file: str
    hits: tuple[Hit, ...]
