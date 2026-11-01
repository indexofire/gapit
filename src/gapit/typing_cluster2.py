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

from typing import Literal

from gapit.cluster import ClusterReport
from gapit.typing_decide import compose_call, decide_scheme
from gapit.typing_engine import gene_calls, score_rules
from gapit.typing_models import TypingDocument
from gapit.typing_results import LocusEvidence, PhenotypeDetail, SchemeCall


def _fallback_locus(report: ClusterReport, scheme_name: str) -> LocusEvidence | None:
    """The best locus in a scheme's family for a fallback call.

    O/K schemes are family-scoped: an O fallback must not report a K locus
    (and vice versa). The report's loci are rank-sorted, so the first member
    of the family is the best evidence even when its coverage is below the
    best-call floor.
    """
    prefix = "OL" if scheme_name == "o_group" else "KL" if scheme_name == "k_group" else None
    if prefix is None:
        return None
    for locus in report.loci:
        if locus.locus.startswith(prefix):
            return LocusEvidence(
                locus=locus.locus,
                coverage_pct=locus.coverage_pct,
                identity_pct=locus.identity_pct,
                genes_present=sum(1 for gene in locus.genes if gene.verdict == "present"),
                genes_partial=sum(1 for gene in locus.genes if gene.verdict == "partial"),
                missing=locus.missing,
            )
    return None


def _scheme_call(
    phenotype: str | None,
    detail: PhenotypeDetail,
    *,
    fallback_locus: LocusEvidence | None = None,
    locus_status: Literal["known", "unknown_locus", "none"] | None = None,
    locus_confidence: Literal["typeable", "untypeable"] | None = None,
) -> SchemeCall:
    """Flatten one decision into a SchemeCall, surfacing the winning
    locus's evidence as the additive ``winner_*`` fields. A profile may
    supply the best family locus for a fallback call, so ``OUT``/``KUT``
    still carries the detected locus and its coverage."""
    winner = detail.winner or fallback_locus
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
        locus_status=locus_status,
        locus_confidence=locus_confidence,
    )


def evaluate_cluster2(
    report: ClusterReport,
    document: TypingDocument,
    *,
    fallback_evidence: bool = False,
    typeable_coverage: float | None = None,
    typeable_identity: float | None = None,
) -> dict[str, SchemeCall]:
    """One SchemeCall per scheme name over the whole loci report (compose
    schemes evaluate after their ingredients, in document order) — the
    evaluation core the gapit.cluster/2 renderers consume. With
    ``fallback_evidence``, a low-confidence fallback still carries the best
    same-family locus so KUT/OUT is distinguishable from "no locus";
    ``typeable_coverage`` / ``typeable_identity`` additionally split the
    evidence into Kaptive-style typeable/untypeable calls."""
    genes = gene_calls(report)
    evaluated: dict[str, SchemeCall] = {}
    for scheme in document.schemes:
        if scheme.compose is not None:
            continue
        phenotype, detail = decide_scheme(scheme, score_rules(report, scheme), genes)
        fallback = (
            _fallback_locus(report, scheme.name)
            if fallback_evidence and detail.confidence == "low" and detail.winner is None
            else None
        )
        winner = (detail.winner or fallback) if fallback_evidence else None
        status: Literal["known", "unknown_locus", "none"] | None = None
        if fallback_evidence:
            if winner is None:
                status = "none"
            elif phenotype is not None and phenotype != scheme.fallback:
                status = "known"
            else:
                status = "unknown_locus"
        confidence: Literal["typeable", "untypeable"] | None = None
        if fallback_evidence and typeable_coverage is not None and typeable_identity is not None:
            confidence = (
                "typeable"
                if winner is not None
                and winner.coverage_pct >= typeable_coverage
                and winner.identity_pct >= typeable_identity
                else "untypeable"
            )
        evaluated[scheme.name] = _scheme_call(
            phenotype,
            detail,
            fallback_locus=fallback,
            locus_status=status,
            locus_confidence=confidence,
        )
    return {
        scheme.name: (
            evaluated[scheme.name]
            if scheme.compose is None
            else compose_call(scheme.compose, evaluated)
        )
        for scheme in document.schemes
    }
