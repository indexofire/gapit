"""gapit.typing evaluation on the gene path — the DEC marker-panel pathway.

A typed GENE database (a FASTA input built with ``db build --typing``)
designates samples through the two-stage pipeline: ``gapit screen -o
result.tsv`` detects the genes, then ``gapit typing result.tsv`` (the
use-case in gapit.typing_ops) folds each FILE's rows into one present
GeneCall per gene name and scores them here — for every scheme, each
weighted_genes/exact_set rule scores over that gene presence, and the
decision layer produces one :class:`SchemeCall` per scheme name, rendered
as gapit.typing_result/1. The cluster path keeps its integrated typing
(gapit.cluster/1's best.phenotype): a locus db has no screen TSV to
re-read, so its designation rides the screen itself.

Stage-2 primitives on this path: ``exact_set`` rules (Doumith/Shigella
marker tables), the scheme-level ``control_gene`` gate and
``unique_group`` mixed-infection override (shared decision layer), rule
``notes``, and ``compose`` schemes rendered from the sibling schemes'
calls after they evaluate.

Feature extraction: a gene counts as present for a rule when its best row
(a hit already cleared the run's minid/mincov filters) clears the rule's
identity_floor — the cluster path's floor semantics, reused verbatim.
Unknown gene references fail at BUILD time (TYPING_UNKNOWN_GENE, the
cluster precedent); gene-kind documents may hold weighted_genes and
exact_set rules only — cluster_match and learned_linear need loci
(TYPING_MALFORMED).
"""

from gapit.cluster import GeneCall
from gapit.errors import DatabaseError
from gapit.typing_decide import decide_scheme
from gapit.typing_engine import score_exact_set, score_weighted_genes
from gapit.typing_models import TypingDocument, TypingScheme, scheme_genes, template_placeholders
from gapit.typing_results import (
    SchemeCall,
    ScoreComponent,
    ScoredRule,
)
from gapit.typing_rules import (
    ClusterMatchRule,
    ExactSetRule,
    LearnedLinearRule,
    WeightedGenesRule,
)

_DETAIL_PLACES = 4


def _unsupported_rule(rule: ClusterMatchRule | LearnedLinearRule) -> DatabaseError:
    return DatabaseError(
        f"typing rule {rule.phenotype!r} ({rule.model}) cannot evaluate against a gene database"
        " (gene databases support weighted_genes and exact_set rules only)",
        code="TYPING_MALFORMED",
        context={"rule": rule.phenotype, "model": rule.model},
    )


def validate_gene_typing(document: TypingDocument, genes: frozenset[str]) -> None:
    """``db build --typing`` on a gene (FASTA) input: every rule must be
    weighted_genes or exact_set, and every referenced gene — rule bodies,
    control genes, unique-group members — must exist in the FASTA records
    (TYPING_UNKNOWN_GENE, the cluster build precedent) — enforced before
    any artifact is written."""
    for scheme in document.schemes:
        for gene_id in scheme_genes(scheme):
            if gene_id not in genes:
                raise DatabaseError(
                    f"typing scheme {scheme.name!r} references unknown gene {gene_id!r}"
                    " (not in the FASTA records)",
                    code="TYPING_UNKNOWN_GENE",
                    context={"gene": gene_id, "scheme": scheme.name, "kind": "gene"},
                )
        for rule in scheme.rules:
            match rule:
                case WeightedGenesRule():
                    names = (*rule.weights, *rule.negative, *rule.require_any)
                case ExactSetRule():
                    names = (*rule.requires, *rule.requires_any, *rule.excludes)
                case ClusterMatchRule() | LearnedLinearRule():
                    raise _unsupported_rule(rule)
            for name in names:
                if name not in genes:
                    raise DatabaseError(
                        f"typing rule {rule.phenotype!r} references unknown gene"
                        f" {name!r} (not in the FASTA records)",
                        code="TYPING_UNKNOWN_GENE",
                        context={"gene": name, "rule": rule.phenotype, "kind": "gene"},
                    )


def score_gene_scheme(scheme: TypingScheme, calls: dict[str, GeneCall]) -> list[ScoredRule]:
    """Every rule of one scheme scored against the gene calls (the gene-path
    twin of typing_engine.score_rules; cluster_match/learned_linear rules
    cannot score here and fail TYPING_MALFORMED)."""
    scored: list[ScoredRule] = []
    for rule in scheme.rules:
        exact = False
        match rule:
            case WeightedGenesRule():
                score, components = score_weighted_genes(rule, calls)
            case ExactSetRule():
                score, components = score_exact_set(rule, calls)
                exact = True
            case ClusterMatchRule() | LearnedLinearRule():
                raise _unsupported_rule(rule)
        scored.append(ScoredRule(rule.phenotype, score, components, tuple(rule.notes), exact))
    return scored


def compose_call(template: str, calls: dict[str, SchemeCall]) -> SchemeCall:
    """Render a compose scheme from its sibling schemes' calls (pure):
    every ``{name}`` placeholder substitutes that scheme's called
    phenotype — fallback strings flow through (a low-confidence ingredient
    makes the composed call low too) — while the FIRST ambiguous
    ingredient (template order) yields a null phenotype carrying that
    ingredient's ambiguous variants and a note naming it. Score is the
    minimum ingredient score; one component per ingredient."""
    names = template_placeholders(template)
    score = round(min(calls[name].score for name in names), _DETAIL_PLACES)
    components = tuple(
        ScoreComponent(name=name, score=round(calls[name].score, _DETAIL_PLACES)) for name in names
    )
    values: dict[str, str] = {}
    for name in names:
        ingredient = calls[name]
        if ingredient.phenotype is None:
            return SchemeCall(
                phenotype=None,
                score=score,
                confidence="ambiguous",
                components=components,
                ambiguous=ingredient.ambiguous,
                notes=(f"ingredient scheme {name!r} is ambiguous",),
            )
        values[name] = ingredient.phenotype
    confidence = "low" if any(calls[name].confidence == "low" for name in names) else "high"
    return SchemeCall(
        phenotype=template.format_map(values),
        score=score,
        confidence=confidence,
        components=components,
    )


def evaluate_gene_calls(
    calls: dict[str, GeneCall], document: TypingDocument
) -> dict[str, SchemeCall]:
    """One SchemeCall per scheme name for one file's gene calls (compose
    schemes evaluate after their ingredients, in document order) — the
    evaluation core ``gapit typing`` runs over each FILE's folded rows."""
    evaluated: dict[str, SchemeCall] = {}
    for scheme in document.schemes:
        if scheme.compose is not None:
            continue
        phenotype, detail = decide_scheme(scheme, score_gene_scheme(scheme, calls), calls)
        evaluated[scheme.name] = SchemeCall(
            phenotype=phenotype,
            score=detail.score,
            confidence=detail.confidence,
            components=detail.components,
            runner_up=detail.runner_up,
            ambiguous=detail.ambiguous,
            notes=detail.notes,
        )
    return {
        scheme.name: (
            evaluated[scheme.name]
            if scheme.compose is None
            else compose_call(scheme.compose, evaluated)
        )
        for scheme in document.schemes
    }
