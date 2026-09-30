"""gapit.typing/1 — the declarative scoring schema of a cluster database.

A typing document is a DATABASE-SIDE artifact (``typing.json`` inside a
cluster db directory, installed via ``gapit db build --typing``): it declares
how hits against a cluster database map to phenotype calls. This module owns
the schema, its validation, and the evaluation RESULT models shared with the
``best`` block of gapit.cluster/1 (stage 3); the arithmetic itself lives in
gapit.typing_engine.

Three rule model kinds exist, discriminated by ``model``:

- ``weighted_genes`` — gene-presence score: positive/negative gene weight
  maps, an identity floor a gene must clear to count, and an any-of gene set
  the locus must match at least once.
- ``cluster_match`` — locus-level component score: weighted coverage,
  identity, and key-genes components, each with its own floor. Component
  weights are documented as summing near 1 but are NOT forced to (the total
  may exceed 1; the cutoff compares the raw score).
- ``learned_linear`` — a trained linear model over named features with a
  bias and free-text training provenance.

Validation failures raise typed DatabaseErrors (TYPING_MALFORMED) carrying
the file path; a missing file is the INPUT_NOT_FOUND input error (the --tsv
precedent).
"""

from pathlib import Path
from typing import Annotated, Literal, NamedTuple

from pydantic import BaseModel, ConfigDict, Field, TypeAdapter, ValidationError

from gapit.errors import DatabaseError, InputError
from gapit.gbfeatures import FeaturesDocument


class WeightedGenesRule(BaseModel, frozen=True):
    """``weighted_genes`` — gene-presence scoring."""

    model: Literal["weighted_genes"]
    phenotype: str = Field(min_length=1)
    weights: dict[str, float]
    negative: dict[str, float] = Field(default_factory=dict)
    identity_floor: float = Field(ge=0.0, le=100.0)
    require_any: list[str] = Field(default_factory=list)


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


class LearnedLinearRule(BaseModel, frozen=True):
    """``learned_linear`` — trained linear model over named features."""

    model: Literal["learned_linear"]
    phenotype: str = Field(min_length=1)
    features: list[str]
    weights: dict[str, float]
    bias: float
    trained_on: dict[str, str]


TypingRule = Annotated[
    WeightedGenesRule | ClusterMatchRule | LearnedLinearRule, Field(discriminator="model")
]


class TypingDocument(BaseModel, frozen=True):
    """gapit.typing/1 — the typing spec document of a cluster database."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.typing/1"] = Field(default="gapit.typing/1", alias="schema")
    rules: tuple[TypingRule, ...] = Field(min_length=1)
    cutoff: float = Field(ge=0.0, le=1.0)
    ambiguity_margin: float = Field(ge=0.0)
    fallback: str


_document_adapter: TypeAdapter[TypingDocument] = TypeAdapter(TypingDocument)


class ScoreComponent(BaseModel, frozen=True):
    """One named contribution inside a phenotype call's score breakdown."""

    name: str
    score: float


class PhenotypeScore(BaseModel, frozen=True):
    """A (phenotype, score) pair: the runner-up or one ambiguous candidate."""

    phenotype: str
    score: float


class PhenotypeDetail(BaseModel, frozen=True):
    """The explainable breakdown behind a phenotype call (the
    ``phenotype_detail`` field of gapit.cluster/1's best block; additive,
    present only when the database carries a typing.json)."""

    score: float
    confidence: Literal["high", "ambiguous", "low"]
    components: tuple[ScoreComponent, ...]
    runner_up: PhenotypeScore | None = None
    ambiguous: tuple[PhenotypeScore, ...] = ()


class ScoredRule(NamedTuple):
    """One evaluated rule: its phenotype, raw score, and component
    breakdown (the decision layer's input)."""

    phenotype: str
    score: float
    components: tuple[ScoreComponent, ...]


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


def validate_references(document: TypingDocument, features: FeaturesDocument) -> None:
    """Every gene/locus id a rule references must exist in the database's
    feature table (TYPING_UNKNOWN_GENE, one code for both kinds) — enforced
    at ``db build --typing`` time and again at screen time."""
    gene_ids = {gene.gene_id for locus in features.loci for gene in locus.genes}
    locus_ids = {locus.id for locus in features.loci}
    for rule in document.rules:
        match rule:
            case WeightedGenesRule():
                referenced = (
                    ("gene", name) for name in (*rule.weights, *rule.negative, *rule.require_any)
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


def _reason(exc: ValidationError) -> str:
    """The validation failure rendered as one line (str(exc) is multiline)."""
    return " ".join(str(exc).split())


def read_typing_document(path: Path) -> TypingDocument:
    """Read and validate a gapit.typing/1 file; a missing file is the typed
    INPUT_NOT_FOUND input error, malformed content the TYPING_MALFORMED
    DatabaseError — both carrying the file in context."""
    if not path.is_file():
        raise InputError(
            f"typing file not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        )
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise DatabaseError(
            f"could not read typing document {path}: {exc}",
            code="TYPING_MALFORMED",
            context={"file": str(path)},
        ) from exc
    try:
        return _document_adapter.validate_json(text)
    except ValidationError as exc:
        raise DatabaseError(
            f"malformed typing document {path}: {_reason(exc)}",
            code="TYPING_MALFORMED",
            context={"file": str(path)},
        ) from exc
