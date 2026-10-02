"""The typing decision layer — what a scheme CALLS.

Split from typing_engine at the 250-LOC ceiling (typing/2 stage 2): this
module owns the decision (ranking, cutoff/ambiguity-margin, and the stage-2
scheme gates), while typing_engine owns how each rule scores. Both
evaluation paths funnel through :func:`decide_scheme`: the cluster path
(:func:`gapit.typing_engine.evaluate_typing`) and the gene path
(:func:`gapit.typing_gene.evaluate_gene_calls`).

Decision pipeline, in order:

1. **Control-gene gate** — a scheme with ``control_gene`` set (prs/ipaH
   style) whose gene is absent by verdict semantics outputs its fallback
   with low confidence and the note ``control gene absent``.
2. **Unique-group override** — when more than one gene of a
   ``unique_group`` is present (two-wzx contamination), the scheme outputs
   its ``mixed_phenotype`` with the present pair listed in ``ambiguous``
   and a mixed-infection note.
3. **Ranked decision** — score desc, declared order on ties, then the
   scheme's cutoff + ambiguity margin. One documented exception: a tie at
   1.0 that involves a satisfied exact_set rule does NOT become
   ambiguous-by-margin; declaration order decides (either side of the tie).
"""

from gapit.cluster import GeneCall
from gapit.typing_models import TypingScheme
from gapit.typing_results import PhenotypeDetail, PhenotypeScore, ScoreComponent, ScoredRule

_DETAIL_PLACES = 4


def _verdict_present(call: GeneCall | None) -> bool:
    """The scheme-level presence predicate: verdict semantics only (a call
    whose verdict is present — rule floors do not exist at scheme level)."""
    return call is not None and call.verdict == "present"


def decide_typing(
    scored: list[ScoredRule], scheme: TypingScheme
) -> tuple[str | None, PhenotypeDetail]:
    """Rank evaluated rules (score desc, declared order on ties) and apply
    the scheme's cutoff + ambiguity-margin decision; returns the phenotype
    string (or None when ambiguous) plus the explainable detail of the
    winning rule."""
    ranked = sorted(enumerate(scored), key=lambda pair: (-pair[1].score, pair[0]))
    top = ranked[0][1]
    second = ranked[1][1] if len(ranked) > 1 else None
    exact_tie = (
        second is not None
        and top.score == 1.0
        and second.score == 1.0
        and (top.exact or second.exact)
    )
    separated = second is None or exact_tie or (top.score - second.score) >= scheme.ambiguity_margin
    if top.score >= scheme.cutoff and separated:
        phenotype, confidence = top.phenotype, "high"
    elif top.score >= scheme.cutoff:
        phenotype, confidence = None, "ambiguous"
    else:
        phenotype, confidence = scheme.fallback, "low"
    return phenotype, PhenotypeDetail(
        score=round(top.score, _DETAIL_PLACES),
        confidence=confidence,
        components=top.components,
        runner_up=(
            None
            if second is None
            else PhenotypeScore(
                phenotype=second.phenotype, score=round(second.score, _DETAIL_PLACES)
            )
        ),
        ambiguous=(
            (
                PhenotypeScore(phenotype=top.phenotype, score=round(top.score, _DETAIL_PLACES)),
                PhenotypeScore(
                    phenotype=second.phenotype, score=round(second.score, _DETAIL_PLACES)
                ),
            )
            if confidence == "ambiguous" and second is not None
            else ()
        ),
        notes=top.notes,
    )


def decide_scheme(
    scheme: TypingScheme, scored: list[ScoredRule], genes: dict[str, GeneCall]
) -> tuple[str | None, PhenotypeDetail]:
    """The scheme-level decision both paths share: control-gene gate first,
    then the unique-group mixed-infection override, then the ranked
    cutoff/margin decision."""
    if scheme.control_gene is not None and not _verdict_present(genes.get(scheme.control_gene)):
        return scheme.fallback, PhenotypeDetail(
            score=0.0, confidence="low", components=(), notes=("control gene absent",)
        )
    for group, members in scheme.unique_group.items():
        present = [member for member in members if _verdict_present(genes.get(member))]
        if len(present) > 1 and scheme.mixed_phenotype is not None:
            return scheme.mixed_phenotype, PhenotypeDetail(
                score=1.0,
                confidence="low",
                components=tuple(ScoreComponent(name=member, score=1.0) for member in present),
                ambiguous=tuple(PhenotypeScore(phenotype=member, score=1.0) for member in present),
                notes=(f"mixed infection in unique group {group!r}",),
            )
    return decide_typing(scored, scheme)
