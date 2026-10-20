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


class LocusEvidence(BaseModel, frozen=True):
    """The winning locus behind a cluster_match call: which locus achieved
    the rule's best component sum, with that locus's own coverage/identity
    and gene-verdict counts. Surfaced per scheme by the gapit.cluster/2
    evaluation; None on every other call path (gate/fallback decisions,
    non-locus rules, the gene path)."""

    locus: str
    coverage_pct: float
    identity_pct: float
    genes_present: int
    genes_partial: int
    missing: tuple[str, ...] = ()


class PhenotypeDetail(BaseModel, frozen=True):
    """The explainable breakdown behind a phenotype call (the
    ``phenotype_detail`` field of gapit.cluster/1's best block; additive,
    present only when the database carries a typing.json).

    ``notes`` (typing/2 stage 2) carries the winning rule's informational
    strings verbatim, or the scheme-gate notes (control gene absent, mixed
    infection); it serializes only when non-empty so pre-stage-2 documents
    keep their byte-identical output. ``winner`` carries the winning
    cluster_match rule's locus evidence (set only on rule-based calls) —
    deliberately NEVER serialized: it is in-memory input for the
    gapit.cluster/2 SchemeCall flattening, so the /1 output stays
    byte-frozen."""

    score: float
    confidence: Literal["high", "ambiguous", "low"]
    components: tuple[ScoreComponent, ...]
    runner_up: PhenotypeScore | None = None
    ambiguous: tuple[PhenotypeScore, ...] = ()
    notes: tuple[str, ...] = ()
    winner: LocusEvidence | None = None

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
    scheme name (the ``gapit typing`` command over a typed gene database),
    and of the gapit.cluster/2 phenotypes map (multi-scheme cluster dbs).

    ``phenotype`` is null on an ambiguous call; the plain serializer keeps
    that null explicit even under a renderer's ``exclude_none``
    (a dropped key would erase the no-call signal), while the optional
    ``runner_up``/``ambiguous``/``notes`` members serialize only when
    set. The ``winner_*`` members (additive, cluster/2 only) flatten the
    winning cluster_match locus's evidence — set together or not at all
    (None on gate/fallback calls, compose schemes, and the whole gene
    path) and serialized only when ``winner_locus`` is set."""

    phenotype: str | None
    score: float
    confidence: Literal["high", "ambiguous", "low"]
    components: tuple[ScoreComponent, ...] = ()
    runner_up: PhenotypeScore | None = None
    ambiguous: tuple[PhenotypeScore, ...] = ()
    notes: tuple[str, ...] = ()
    winner_locus: str | None = None
    winner_coverage_pct: float | None = None
    winner_identity_pct: float | None = None
    winner_present: int | None = None
    winner_partial: int | None = None
    winner_missing: tuple[str, ...] | None = None

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
        if self.winner_locus is not None:
            data["winner_locus"] = self.winner_locus
            data["winner_coverage_pct"] = self.winner_coverage_pct
            data["winner_identity_pct"] = self.winner_identity_pct
            data["winner_present"] = self.winner_present
            data["winner_partial"] = self.winner_partial
            data["winner_missing"] = list(self.winner_missing or ())
        return data


class ScoredRule(NamedTuple):
    """One evaluated rule: its phenotype, raw score, component breakdown
    (the decision layer's input), and — stage 2 — the rule's verbatim
    ``notes`` plus ``exact`` (True only for exact_set rules, the flag
    behind the satisfied-1.0 tie exception). ``winner`` (cluster/2) is the
    best locus's evidence a cluster_match rule scored against, None for
    every other rule kind."""

    phenotype: str
    score: float
    components: tuple[ScoreComponent, ...]
    notes: tuple[str, ...] = ()
    exact: bool = False
    winner: LocusEvidence | None = None
