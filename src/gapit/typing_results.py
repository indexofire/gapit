"""gapit.typing evaluation RESULT models (the output side of typing).

Shared by both evaluation paths: PhenotypeDetail annotates gapit.cluster/1's
``best`` block (typed cluster dbs), SchemeCall populates the
``files[].phenotypes`` object of gapit.typing_result/1 (the ``gapit
typing`` command over typed gene dbs), and
ScoreComponent/PhenotypeScore/ScoredRule are the explainable-breakdown
building blocks the decision layer in typing_engine/typing_gene consumes.
Pure pydantic value models — no I/O, no arithmetic.
"""

from typing import Literal, NamedTuple

from pydantic import BaseModel, model_serializer


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
    present only when the database carries a typing.json).

    ``notes`` (typing/2 stage 2) carries the winning rule's informational
    strings verbatim, or the scheme-gate notes (control gene absent, mixed
    infection); it serializes only when non-empty so pre-stage-2 documents
    keep their byte-identical output."""

    score: float
    confidence: Literal["high", "ambiguous", "low"]
    components: tuple[ScoreComponent, ...]
    runner_up: PhenotypeScore | None = None
    ambiguous: tuple[PhenotypeScore, ...] = ()
    notes: tuple[str, ...] = ()

    @model_serializer
    def _serialize(self) -> dict[str, object]:
        data: dict[str, object] = {
            "score": self.score,
            "confidence": self.confidence,
            "components": [component.model_dump(mode="json") for component in self.components],
            "runner_up": (
                None if self.runner_up is None else self.runner_up.model_dump(mode="json")
            ),
            "ambiguous": [entry.model_dump(mode="json") for entry in self.ambiguous],
        }
        if self.notes:
            data["notes"] = list(self.notes)
        return data


class SchemeCall(BaseModel, frozen=True):
    """One scheme's call for a screened file — one entry of the
    ``files[].phenotypes`` object of gapit.typing_result/1, keyed by
    scheme name (the ``gapit typing`` command over a typed gene database).

    ``phenotype`` is null on an ambiguous call; the plain serializer keeps
    that null explicit even under a renderer's ``exclude_none``
    (a dropped key would erase the no-call signal), while the optional
    ``runner_up``/``ambiguous``/``notes`` members serialize only when
    set."""

    phenotype: str | None
    score: float
    confidence: Literal["high", "ambiguous", "low"]
    components: tuple[ScoreComponent, ...] = ()
    runner_up: PhenotypeScore | None = None
    ambiguous: tuple[PhenotypeScore, ...] = ()
    notes: tuple[str, ...] = ()

    @model_serializer
    def _serialize(self) -> dict[str, object]:
        data: dict[str, object] = {
            "phenotype": self.phenotype,
            "score": self.score,
            "confidence": self.confidence,
            "components": [component.model_dump(mode="json") for component in self.components],
        }
        if self.runner_up is not None:
            data["runner_up"] = self.runner_up.model_dump(mode="json")
        if self.ambiguous:
            data["ambiguous"] = [entry.model_dump(mode="json") for entry in self.ambiguous]
        if self.notes:
            data["notes"] = list(self.notes)
        return data


class ScoredRule(NamedTuple):
    """One evaluated rule: its phenotype, raw score, component breakdown
    (the decision layer's input), and — stage 2 — the rule's verbatim
    ``notes`` plus ``exact`` (True only for exact_set rules, the flag
    behind the satisfied-1.0 tie exception)."""

    phenotype: str
    score: float
    components: tuple[ScoreComponent, ...]
    notes: tuple[str, ...] = ()
    exact: bool = False
