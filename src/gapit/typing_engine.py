"""The gapit.typing evaluation engine — the cluster path (stage 3).

Pure functions: one typing SCHEME (:class:`gapit.typing_models.TypingScheme`
— a gapit.typing/2 scheme, or the single scheme a gapit.typing/1 document
degrades to) plus one file's cluster results (stage 2's ClusterReport)
become a phenotype call with an explainable score breakdown. The three rule
models score independently —
``weighted_genes`` normalizes gene-presence weights to [0, 1],
``cluster_match`` sums raw weighted locus components (per covered locus, best
locus wins, rank-order tie-break; a component below its floor zeroes the
locus's whole score), ``learned_linear`` squashes weighted named features
through a sigmoid — and the decision layer ranks rules by score desc
(declared order on ties): top >= cutoff AND (top - second) >=
ambiguity_margin calls the phenotype with ``high`` confidence; top >= cutoff
with a thinner separation is ``ambiguous`` (phenotype null, the top two
listed); otherwise the scheme's fallback with ``low`` confidence.

The gene-kind pathway (gapit.typing_result/1 calls over weighted_genes
rules) lives in gapit.typing_gene and reuses this module's scorers verbatim.
Reference validation (TYPING_UNKNOWN_GENE) lives in typing_models beside the
schema; loading (typing.json beside sequences) in gapit.cluster beside
load_features. Reports without a best call pass through unannotated.
"""

import math

from gapit.cluster import ClusterReport, GeneCall, LocusCall
from gapit.errors import DatabaseError
from gapit.typing_decide import decide_scheme
from gapit.typing_models import TypingScheme
from gapit.typing_results import ScoreComponent, ScoredRule
from gapit.typing_rules import (
    ClusterMatchRule,
    ExactSetRule,
    LearnedLinearRule,
    WeightedGenesRule,
    parse_feature_name,
)

_DETAIL_PLACES = 4


def _clamp01(score: float) -> float:
    return min(1.0, max(0.0, score))


def _component(name: str, score: float) -> ScoreComponent:
    return ScoreComponent(name=name, score=round(score, _DETAIL_PLACES))


def score_weighted_genes(
    rule: WeightedGenesRule, genes: dict[str, GeneCall]
) -> tuple[float, tuple[ScoreComponent, ...]]:
    """Gene-presence score normalized to [0, 1] by the sum of positive
    weights (raw sum when that is 0). A gene counts when its verdict is
    present, identity clears the rule's identity floor, and — when the
    optional coverage_floor is set — coverage clears that too; an
    unsatisfied require_any zeroes the rule outright."""

    def present(gene_id: str) -> bool:
        call = genes.get(gene_id)
        return (
            call is not None
            and call.verdict == "present"
            and call.identity_pct >= rule.identity_floor
            and (rule.coverage_floor is None or call.coverage_pct >= rule.coverage_floor)
        )

    if rule.require_any and not any(present(gene_id) for gene_id in rule.require_any):
        return 0.0, ()
    denominator = sum(rule.weights.values())
    contributions = [
        (gene_id, weight if present(gene_id) else 0.0)
        for gene_id, weight in (*rule.weights.items(), *rule.negative.items())
    ]
    raw = sum(score for _, score in contributions)
    scaled = raw / denominator if denominator > 0 else raw
    return _clamp01(scaled), tuple(
        _component(gene_id, score / denominator if denominator > 0 else score)
        for gene_id, score in contributions
    )


def score_exact_set(
    rule: ExactSetRule, genes: dict[str, GeneCall]
) -> tuple[float, tuple[ScoreComponent, ...]]:
    """Deterministic boolean gene-set match (Doumith/Shigella-style tables):
    1.0 iff every ``requires`` gene is present, at least one
    ``requires_any`` gene is present (when that set is non-empty), AND
    every ``excludes`` gene is absent, else 0.0. Presence is the
    weighted_genes floor semantics with both floors defaulting to 90; each
    component carries that gene's constraint satisfaction (1.0 met / 0.0
    missed), required genes first, then the any-of set, then excludes."""

    def satisfied(gene_id: str) -> bool:
        call = genes.get(gene_id)
        return (
            call is not None
            and call.verdict == "present"
            and call.identity_pct >= rule.identity_floor
            and call.coverage_pct >= rule.coverage_floor
        )

    ok = (
        all(satisfied(gene_id) for gene_id in rule.requires)
        and (not rule.requires_any or any(satisfied(gene_id) for gene_id in rule.requires_any))
        and all(not satisfied(gene_id) for gene_id in rule.excludes)
    )
    components = (
        tuple(_component(gene_id, 1.0 if satisfied(gene_id) else 0.0) for gene_id in rule.requires)
        + tuple(
            _component(gene_id, 1.0 if satisfied(gene_id) else 0.0) for gene_id in rule.requires_any
        )
        + tuple(
            _component(gene_id, 1.0 if not satisfied(gene_id) else 0.0) for gene_id in rule.excludes
        )
    )
    return (1.0 if ok else 0.0), components


def score_cluster_match(
    rule: ClusterMatchRule, report: ClusterReport
) -> tuple[float, tuple[ScoreComponent, ...]]:
    """Best covered locus's weighted component sum (raw, not normalized):
    below-floor coverage or identity zeroes THAT locus's score; key genes
    are the fraction of the listed ids present inside the locus (an empty
    list scores 0). Ties keep the earlier rank."""
    winner: tuple[float, tuple[ScoreComponent, ...]] | None = None
    for locus in report.loci:
        if locus.coverage_pct < rule.coverage.floor or locus.identity_pct < rule.identity.floor:
            continue
        calls = {call.gene_id: call for call in locus.genes}
        key_fraction = (
            sum(
                1
                for gene_id in rule.key_genes.genes
                if calls.get(gene_id) is not None and calls[gene_id].verdict == "present"
            )
            / len(rule.key_genes.genes)
            if rule.key_genes.genes
            else 0.0
        )
        coverage = locus.coverage_pct / 100.0
        identity = locus.identity_pct / 100.0
        components = (
            _component("coverage", rule.coverage.weight * coverage),
            _component("identity", rule.identity.weight * identity),
            _component("key_genes", rule.key_genes.weight * key_fraction),
        )
        total = sum(component.score for component in components)
        if winner is None or total > winner[0]:
            winner = (total, components)
    return winner if winner is not None else (0.0, ())


def score_learned_linear(
    rule: LearnedLinearRule, report: ClusterReport, genes: dict[str, GeneCall]
) -> tuple[float, tuple[ScoreComponent, ...]]:
    """sigmoid(sum(weight_i * feature_i) + bias); every listed feature must
    have a weight (TYPING_MALFORMED otherwise) and reads fractions — gene
    cov/ident/present (present = verdict present), cluster coverage/identity
    (0 when the locus was never covered)."""

    def value(kind: str, name: str, metric: str) -> float:
        if kind == "gene":
            call = genes.get(name)
            if call is None:
                return 0.0
            match metric:
                case "cov":
                    return call.coverage_pct / 100.0
                case "ident":
                    return call.identity_pct / 100.0
                case "present":
                    return 1.0 if call.verdict == "present" else 0.0
                case _:
                    raise _malformed_metric(kind, name, metric)
        locus: LocusCall | None = next(
            (entry for entry in report.loci if entry.locus == name), None
        )
        if locus is None:
            return 0.0
        match metric:
            case "coverage":
                return locus.coverage_pct / 100.0
            case "identity":
                return locus.identity_pct / 100.0
            case _:
                raise _malformed_metric(kind, name, metric)

    total = rule.bias
    components = [_component("bias", rule.bias)]
    for feature in rule.features:
        weight = rule.weights.get(feature)
        if weight is None:
            raise DatabaseError(
                f"learned_linear feature {feature!r} has no weight",
                code="TYPING_MALFORMED",
                context={"feature": feature},
            )
        kind, name, metric = parse_feature_name(feature)
        contribution = weight * value(kind, name, metric)
        components.append(_component(feature, contribution))
        total += contribution
    if total >= 0:
        return 1.0 / (1.0 + math.exp(-total)), tuple(components)
    exp_total = math.exp(total)
    return exp_total / (1.0 + exp_total), tuple(components)


def _malformed_metric(kind: str, name: str, metric: str) -> DatabaseError:
    return DatabaseError(
        f"unknown learned_linear feature metric {metric!r} in {kind}:{name}:{metric}",
        code="TYPING_MALFORMED",
        context={"feature": f"{kind}:{name}:{metric}"},
    )


def gene_calls(report: ClusterReport) -> dict[str, GeneCall]:
    """One call per gene id across all loci (the first locus wins — the
    fold score_rules scores against)."""
    genes: dict[str, GeneCall] = {}
    for locus in report.loci:
        for call in locus.genes:
            genes.setdefault(call.gene_id, call)
    return genes


def score_rules(report: ClusterReport, scheme: TypingScheme) -> list[ScoredRule]:
    """Every rule of one scheme scored against one report (the shared
    dispatch behind evaluate_typing; the calibration harness consumes it to
    read the score of ANY phenotype's rule, not just the winner)."""
    genes = gene_calls(report)
    scored: list[ScoredRule] = []
    for rule in scheme.rules:
        exact = False
        match rule:
            case WeightedGenesRule():
                score, components = score_weighted_genes(rule, genes)
            case ExactSetRule():
                score, components = score_exact_set(rule, genes)
                exact = True
            case ClusterMatchRule():
                score, components = score_cluster_match(rule, report)
            case LearnedLinearRule():
                score, components = score_learned_linear(rule, report, genes)
        scored.append(ScoredRule(rule.phenotype, score, components, tuple(rule.notes), exact))
    return scored


def evaluate_typing(report: ClusterReport, scheme: TypingScheme) -> ClusterReport:
    """Annotate a file's report: best gains phenotype + phenotype_detail
    (reports without a best call pass through untouched)."""
    if report.best is None:
        return report
    phenotype, detail = decide_scheme(scheme, score_rules(report, scheme), gene_calls(report))
    return report.model_copy(
        update={
            "best": report.best.model_copy(
                update={"phenotype": phenotype, "phenotype_detail": detail}
            )
        }
    )
