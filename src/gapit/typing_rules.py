"""gapit.typing rule models — the per-phenotype scoring units.

Split from typing_models.py at the 250-LOC ceiling (typing/2 stage 2): this
module owns the RULE layer only — the four ``model`` kinds, their optional
floors and ``notes`` channels, the discriminated :data:`TypingRule` union,
and the rule-level reference checks. The document container, its parsing,
and the scheme-level fields stay in typing_models.

Rule kinds (discriminated by ``model``):

- ``weighted_genes`` — gene-presence score: positive/negative gene weight
  maps, an identity floor a gene must clear to count, an any-of gene set
  the locus must match at least once, and (stage 2) an optional
  ``coverage_floor`` that additionally gates a gene's coverage.
- ``exact_set`` (stage 2) — deterministic boolean gene-set match
  (Doumith/Shigella-style marker tables): 1.0 iff every ``requires`` gene
  is present AND at least one ``requires_any`` gene is present (when that
  set is non-empty) AND every ``excludes`` gene is absent, else 0.0. The
  three gene sets must be pairwise disjoint and not all empty. Presence
  uses the weighted_genes floor semantics with both floors defaulting to
  90. A satisfied exact_set at 1.0 tied with another rule at 1.0 does NOT
  become ambiguous-by-margin: declaration order decides (the decision
  layer's documented exception).
- ``cluster_match`` — locus-level component score (cluster dbs only).
- ``learned_linear`` — a trained linear model over named features.

Every rule also carries optional ``notes`` (stage 2): free-text strings
surfaced verbatim on the phenotype call when the rule wins
(pINV-ipaB/phase-OFF-style informational flags).
"""

from collections.abc import Iterator
from typing import Annotated, Literal

from pydantic import BaseModel, Field, model_validator

from gapit.errors import DatabaseError


class WeightedGenesRule(BaseModel, frozen=True):
    """``weighted_genes`` — gene-presence scoring. ``coverage_floor``
    (optional) additionally requires a gene's coverage to clear the floor
    before it counts; unset means no gene-coverage gate beyond the engine's
    verdict semantics."""

    model: Literal["weighted_genes"]
    phenotype: str = Field(min_length=1)
    weights: dict[str, float]
    negative: dict[str, float] = Field(default_factory=dict)
    identity_floor: float = Field(ge=0.0, le=100.0)
    require_any: list[str] = Field(default_factory=list)
    coverage_floor: float | None = Field(default=None, ge=0.0, le=100.0)
    notes: list[str] = Field(default_factory=list)


class ExactSetRule(BaseModel, frozen=True):
    """``exact_set`` — deterministic boolean gene-set match: score 1.0 iff
    every ``requires`` gene is present AND at least one ``requires_any``
    gene is present (when that set is non-empty) AND every ``excludes``
    gene is absent, else 0.0. Presence gates on both floors (90/90 by
    default, the engine's verdict defaults)."""

    model: Literal["exact_set"]
    phenotype: str = Field(min_length=1)
    requires: list[str] = Field(default_factory=list)
    requires_any: list[str] = Field(default_factory=list)
    excludes: list[str] = Field(default_factory=list)
    identity_floor: float = Field(default=90.0, ge=0.0, le=100.0)
    coverage_floor: float = Field(default=90.0, ge=0.0, le=100.0)
    notes: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def _disjoint(self) -> "ExactSetRule":
        if not (self.requires or self.requires_any or self.excludes):
            raise ValueError(
                "exact_set rule declares no genes (requires, requires_any, and excludes"
                " are all empty)"
            )
        overlaps = (
            ("required and in requires_any", set(self.requires) & set(self.requires_any)),
            ("required and excluded", set(self.requires) & set(self.excludes)),
            ("in requires_any and excluded", set(self.requires_any) & set(self.excludes)),
        )
        for description, overlap in overlaps:
            if overlap:
                raise ValueError(f"exact_set genes are both {description}: {sorted(overlap)}")
        return self


class CoverageComponent(BaseModel, frozen=True):
    weight: float
    floor: float = Field(ge=0.0, le=100.0)


class IdentityComponent(BaseModel, frozen=True):
    weight: float
    floor: float = Field(ge=0.0, le=100.0)


class KeyGenesComponent(BaseModel, frozen=True):
    weight: float
    genes: list[str]


class ClusterMatchRule(BaseModel, frozen=True):
    """``cluster_match`` — locus-level component scoring."""

    model: Literal["cluster_match"]
    phenotype: str = Field(min_length=1)
    coverage: CoverageComponent
    identity: IdentityComponent
    key_genes: KeyGenesComponent
    notes: list[str] = Field(default_factory=list)


class LearnedLinearRule(BaseModel, frozen=True):
    """``learned_linear`` — trained linear model over named features."""

    model: Literal["learned_linear"]
    phenotype: str = Field(min_length=1)
    features: list[str]
    weights: dict[str, float]
    bias: float
    trained_on: dict[str, str]
    notes: list[str] = Field(default_factory=list)


TypingRule = Annotated[
    WeightedGenesRule | ExactSetRule | ClusterMatchRule | LearnedLinearRule,
    Field(discriminator="model"),
]


def parse_feature_name(feature: str) -> tuple[str, str, str]:
    """``gene:<id>:<cov|ident|present>`` / ``cluster:<locus>:<coverage|identity>``
    -> (kind, id, metric); anything else is TYPING_MALFORMED."""
    parts = feature.split(":")
    if len(parts) != 3 or not all(parts) or parts[0] not in ("gene", "cluster"):
        raise DatabaseError(
            f"malformed learned_linear feature {feature!r}"
            " (expected gene:<id>:<cov|ident|present> or cluster:<locus>:<coverage|identity>)",
            code="TYPING_MALFORMED",
            context={"feature": feature},
        )
    return parts[0], parts[1], parts[2]


def check_rule_references(
    rule: TypingRule, gene_ids: frozenset[str] | set[str], locus_ids: frozenset[str] | set[str]
) -> None:
    """Every gene/locus id a rule references must exist in the database's
    feature table (TYPING_UNKNOWN_GENE, one code for both kinds)."""
    referenced: Iterator[tuple[str, str]]
    match rule:
        case WeightedGenesRule():
            referenced = (
                ("gene", name) for name in (*rule.weights, *rule.negative, *rule.require_any)
            )
        case ExactSetRule():
            referenced = (
                ("gene", name) for name in (*rule.requires, *rule.requires_any, *rule.excludes)
            )
        case ClusterMatchRule():
            referenced = (("gene", name) for name in rule.key_genes.genes)
        case LearnedLinearRule():
            referenced = (
                (parse_feature_name(feature)[0], parse_feature_name(feature)[1])
                for feature in rule.features
            )
    for kind, name in referenced:
        known = gene_ids if kind == "gene" else locus_ids
        if name not in known:
            key = "gene" if kind == "gene" else "locus"
            raise DatabaseError(
                f"typing rule {rule.phenotype!r} references unknown {kind} {name!r}"
                " (not in the database features)",
                code="TYPING_UNKNOWN_GENE",
                context={key: name, "rule": rule.phenotype, "kind": kind},
            )
