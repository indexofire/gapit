"""Evaluate gapit.serotype/1 rules over one file's gene calls."""

from __future__ import annotations

from collections.abc import Mapping

from gapit.cluster import GeneCall
from gapit.serotype_models import SerotypeDocument, SerotypeRule
from gapit.typing_results import PhenotypeScore, SchemeCall, ScoreComponent


def _present(call: GeneCall | None, coverage: float, identity: float) -> bool:
    return call is not None and call.coverage_pct >= coverage and call.identity_pct >= identity


def _rule_matches(rule: SerotypeRule, calls: Mapping[str, GeneCall]) -> tuple[str, ...] | None:
    """The ordered positive markers when the rule matches, else None."""
    positives: list[str] = []
    for gene in rule.requires:
        if not _present(calls.get(gene), rule.coverage_floor, rule.identity_floor):
            return None
        positives.append(gene)
    if rule.requires_any:
        any_hit = next(
            (
                gene
                for gene in rule.requires_any
                if _present(calls.get(gene), rule.coverage_floor, rule.identity_floor)
            ),
            None,
        )
        if any_hit is None:
            return None
        positives.append(any_hit)
    for gene in rule.excludes:
        if _present(calls.get(gene), rule.coverage_floor, rule.identity_floor):
            return None
    for gene, ceiling in rule.max_coverage.items():
        call = calls.get(gene)
        if call is not None and call.coverage_pct >= ceiling:
            return None
    for gene, minimum in rule.min_alignment_start_pct.items():
        call = calls.get(gene)
        if call is not None and (
            call.alignment_start_pct is None or call.alignment_start_pct < minimum
        ):
            return None
    for gene, maximum in rule.max_alignment_start_pct.items():
        call = calls.get(gene)
        if call is not None and (
            call.alignment_start_pct is None or call.alignment_start_pct > maximum
        ):
            return None
    return tuple(positives)


def _control_present(calls: Mapping[str, GeneCall], document: SerotypeDocument) -> bool:
    if document.control_gene is None:
        return True
    return _present(
        calls.get(document.control_gene),
        document.control_coverage_floor,
        document.control_identity_floor,
    )


def evaluate_serotype(calls: Mapping[str, GeneCall], document: SerotypeDocument) -> SchemeCall:
    """One ordered marker-set decision, surfacing every matching rule.

    A rule whose ``control_exempt`` is set can still fire when the scheme
    control gene is absent (the ShigaTyper S. boydii 13 precedent). Multiple
    distinct phenotypes are an ambiguous call, never declaration-order luck.
    """
    control = _control_present(calls, document)
    matched: list[tuple[SerotypeRule, tuple[str, ...]]] = []
    for rule in document.rules:
        if not control and not rule.control_exempt:
            continue
        positive = _rule_matches(rule, calls)
        if positive is not None:
            matched.append((rule, positive))
    if matched:
        phenotypes = {rule.phenotype for rule, _ in matched}
        positives = sorted({gene for _, genes in matched for gene in genes})
        components = tuple(ScoreComponent(name=f"marker:{gene}", score=1.0) for gene in positives)
        if len(phenotypes) == 1:
            rule, genes = matched[0]
            notes = list(rule.notes)
            notes.append("markers: " + ",".join(genes))
            if not control:
                notes.append("control gene absent; control_exempt rule used")
            return SchemeCall(
                phenotype=rule.phenotype,
                score=1.0,
                confidence="high",
                components=components,
                notes=tuple(notes),
            )
        candidates = tuple(
            PhenotypeScore(phenotype=phenotype, score=1.0) for phenotype in sorted(phenotypes)
        )
        return SchemeCall(
            phenotype=None,
            score=1.0,
            confidence="ambiguous",
            components=components,
            ambiguous=candidates,
            notes=(document.ambiguous, "candidates: " + ",".join(sorted(phenotypes))),
        )
    if not control:
        return SchemeCall(
            phenotype=document.fallback,
            score=0.0,
            confidence="low",
            notes=("control gene absent",),
        )
    return SchemeCall(
        phenotype=document.unresolved,
        score=0.0,
        confidence="low",
        notes=("no marker set matched",),
    )


__all__ = ["evaluate_serotype"]
