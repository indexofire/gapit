"""gapit.typing evaluation on the multi-scheme cluster path (cluster/2).

A cluster database whose typing document holds MORE than one scheme
(:func:`gapit.typing_models.uses_cluster_v2`) renders gapit.cluster/2:
every scheme is called — not just the best-locus annotation the /1 path
makes — and the per-scheme :class:`~gapit.typing_results.SchemeCall` map
rides the report (``ClusterReport.phenotypes``, insertion order = scheme
document order). Evaluation deliberately skips the /1 best-is-None early
return: schemes score over the full loci report (empty loci leave every
rule at 0, so each scheme falls back), and compose schemes render from
their sibling ingredients' calls through the shared decision layer.
"""

from gapit.cluster import ClusterReport
from gapit.typing_decide import compose_call, decide_scheme
from gapit.typing_engine import gene_calls, score_rules
from gapit.typing_models import TypingDocument
from gapit.typing_results import PhenotypeDetail, SchemeCall


def _scheme_call(phenotype: str | None, detail: PhenotypeDetail) -> SchemeCall:
    """Flatten one decision into a SchemeCall, surfacing the winning
    locus's evidence as the additive ``winner_*`` fields (all None when
    the decision is gate- or fallback-based — evidence only exists for a
    rule-based call)."""
    winner = detail.winner
    return SchemeCall(
        phenotype=phenotype,
        score=detail.score,
        confidence=detail.confidence,
        components=detail.components,
        runner_up=detail.runner_up,
        ambiguous=detail.ambiguous,
        notes=detail.notes,
        winner_locus=winner.locus if winner is not None else None,
        winner_coverage_pct=round(winner.coverage_pct, 2) if winner is not None else None,
        winner_identity_pct=round(winner.identity_pct, 2) if winner is not None else None,
        winner_present=winner.genes_present if winner is not None else None,
        winner_partial=winner.genes_partial if winner is not None else None,
        winner_missing=winner.missing if winner is not None else None,
    )


def evaluate_cluster2(report: ClusterReport, document: TypingDocument) -> dict[str, SchemeCall]:
    """One SchemeCall per scheme name over the whole loci report (compose
    schemes evaluate after their ingredients, in document order) — the
    evaluation core the gapit.cluster/2 renderers consume."""
    genes = gene_calls(report)
    evaluated: dict[str, SchemeCall] = {}
    for scheme in document.schemes:
        if scheme.compose is not None:
            continue
        phenotype, detail = decide_scheme(scheme, score_rules(report, scheme), genes)
        evaluated[scheme.name] = _scheme_call(phenotype, detail)
    return {
        scheme.name: (
            evaluated[scheme.name]
            if scheme.compose is None
            else compose_call(scheme.compose, evaluated)
        )
        for scheme in document.schemes
    }
