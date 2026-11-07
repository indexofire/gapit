"""gspec.serotype/1 — marker-set serotype profile documents.

This is deliberately separate from gspec.typing: a serotype call combines
several antigen/POAC/virulence markers and often needs an ordered fallback
and ambiguity policy. The first consumer is the Shigella marker panel
(``shi_sero/serotype.json``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator

from gspec.errors import DatabaseError
from gspec.schema_ids import Canonical

SEROTYPE_FILENAME = "serotype.json"


class SerotypeRule(BaseModel, frozen=True):
    """One ordered marker-set rule; every matched rule is surfaced.

    ``max_coverage`` is a per-marker ceiling (strictly below the value); it
    expresses the S. boydii 6 IS629-interrupted ``wbaM`` case without turning
    the whole rule into a low-coverage call.  ``min_alignment_start_pct`` and
    ``max_alignment_start_pct`` gate a marker's alignment start (0 = the 5'
    end is covered), separating a 5'-interrupted allele from a 3'-truncated
    assembly.  A rule that declares a start gate requires alignment evidence;
    reads tables carry none, so such rules do not match there.
    """

    phenotype: str
    requires: tuple[str, ...] = ()
    requires_any: tuple[str, ...] = ()
    excludes: tuple[str, ...] = ()
    max_coverage: dict[str, float] = Field(
        default_factory=dict,
        description="Per-marker coverage ceiling; the rule fails at or above this percentage.",
    )
    min_alignment_start_pct: dict[str, float] = Field(
        default_factory=dict,
        description=(
            "Per-marker minimum alignment start as percent of marker length"
            " (0 = the 5' end is covered)."
        ),
    )
    max_alignment_start_pct: dict[str, float] = Field(
        default_factory=dict,
        description=(
            "Per-marker maximum alignment start as percent of marker length"
            " (0 = the 5' end is covered)."
        ),
    )
    coverage_floor: float = Field(default=50.0, ge=0.0, le=100.0)
    identity_floor: float = Field(default=80.0, ge=0.0, le=100.0)
    control_exempt: bool = False
    notes: tuple[str, ...] = ()

    @model_validator(mode="after")
    def _has_positive_marker(self) -> SerotypeRule:
        if not (self.requires or self.requires_any):
            raise ValueError("serotype rule needs requires or requires_any")
        if set(self.requires) & set(self.requires_any):
            raise ValueError("requires and requires_any overlap")
        if set(self.requires) & set(self.excludes):
            raise ValueError("requires and excludes overlap")
        if set(self.requires_any) & set(self.excludes):
            raise ValueError("requires_any and excludes overlap")
        positives = {*self.requires, *self.requires_any}
        for name, gates in (
            ("max_coverage", self.max_coverage),
            ("min_alignment_start_pct", self.min_alignment_start_pct),
            ("max_alignment_start_pct", self.max_alignment_start_pct),
        ):
            for gene, value in gates.items():
                if gene not in positives:
                    raise ValueError(f"{name} gene {gene!r} is not a positive marker")
                if not 0.0 <= value <= 100.0:
                    raise ValueError(f"{name}[{gene!r}] must be a percentage, got {value}")
        return self


class SerotypeDocument(BaseModel, frozen=True):
    """gspec.serotype/1 — one database's ordered serotype rules."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Annotated[Literal["gspec.serotype/1"], Canonical] = Field(
        default="gspec.serotype/1", alias="schema"
    )
    name: str
    control_gene: str | None = None
    control_coverage_floor: float = Field(default=50.0, ge=0.0, le=100.0)
    control_identity_floor: float = Field(default=80.0, ge=0.0, le=100.0)
    fallback: str = "Not Shigella/EIEC"
    ambiguous: str = "No prediction (multiple marker sets)"
    unresolved: str = "Shigella spp. (unresolved)"
    rules: tuple[SerotypeRule, ...]


def read_serotype_document(path: Path) -> SerotypeDocument:
    """Read and validate a serotype profile; malformed content is typed."""
    try:
        return SerotypeDocument.model_validate_json(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError) as exc:
        raise DatabaseError(
            f"cannot read serotype profile {path}: {exc}",
            code="SEROTYPE_MALFORMED",
            context={"file": str(path)},
        ) from exc
    except ValidationError as exc:
        raise DatabaseError(
            f"malformed serotype profile {path}: {exc}",
            code="SEROTYPE_MALFORMED",
            context={"file": str(path)},
        ) from exc


def load_serotype(db_dir: Path) -> SerotypeDocument | None:
    """The optional serotype profile next to a gene database."""
    path = db_dir / SEROTYPE_FILENAME
    return read_serotype_document(path) if path.is_file() else None


def validate_serotype_references(document: SerotypeDocument, gene_ids: frozenset[str]) -> None:
    """Every marker a rule references must exist in the database."""
    for rule in document.rules:
        for gene in (
            *rule.requires,
            *rule.requires_any,
            *rule.excludes,
            *rule.max_coverage,
            *rule.min_alignment_start_pct,
            *rule.max_alignment_start_pct,
        ):
            if gene not in gene_ids:
                raise DatabaseError(
                    f"serotype rule {rule.phenotype!r} references unknown marker {gene!r}",
                    code="SEROTYPE_UNKNOWN_GENE",
                    context={"phenotype": rule.phenotype, "gene": gene},
                )


__all__ = [
    "SEROTYPE_FILENAME",
    "SerotypeDocument",
    "SerotypeRule",
    "load_serotype",
    "read_serotype_document",
    "validate_serotype_references",
]
