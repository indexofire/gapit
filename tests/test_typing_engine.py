"""Typing-engine unit suites (stage 3): the pure gapit.typing/1 evaluators.

No I/O anywhere: ClusterReport/LocusCall/GeneCall values are built directly,
rules come from typing_models, and every score/decision is asserted on the
exact arithmetic (the SPEC-style formulas, not implementation echoes). The
CLI/golden integration for typed databases lives in test_cluster_typing.py.
"""

import math
from pathlib import Path

import pytest
from pydantic import TypeAdapter

from gapit.cluster import BestCall, ClusterReport, GeneCall, LocusCall, load_typing
from gapit.cluster_math import Verdict
from gapit.db import Database
from gapit.errors import DatabaseError
from gapit.gbfeatures import FeaturesDocument, GeneFeature, LocusFeatures
from gapit.typing_engine import (
    decide_typing,
    evaluate_typing,
    score_cluster_match,
    score_learned_linear,
    score_weighted_genes,
)
from gapit.typing_models import (
    ClusterMatchRule,
    CoverageComponent,
    IdentityComponent,
    KeyGenesComponent,
    LearnedLinearRule,
    ScoreComponent,
    ScoredRule,
    TypingDocument,
    WeightedGenesRule,
    validate_references,
)


def gene(
    gene_id: str, verdict: Verdict = "present", coverage: float = 100.0, identity: float = 100.0
) -> GeneCall:
    return GeneCall(
        gene_id=gene_id,
        start=1,
        end=100,
        strand="+",
        coverage_pct=coverage,
        identity_pct=identity,
        verdict=verdict,
    )


def locus_call(
    locus_id: str,
    genes: list[GeneCall],
    coverage: float = 100.0,
    identity: float = 100.0,
    rank: int = 1,
) -> LocusCall:
    return LocusCall(
        locus=locus_id,
        label=f"label {locus_id}",
        type=f"KL_{locus_id}",
        coverage_pct=coverage,
        identity_pct=identity,
        rank=rank,
        genes=tuple(genes),
        missing=tuple(g.gene_id for g in genes if g.verdict != "present"),
    )


def report(loci: list[LocusCall]) -> ClusterReport:
    best: BestCall | None = None
    if loci and loci[0].coverage_pct >= 96.0:
        top = loci[0]
        best = BestCall(
            locus=top.locus,
            label=top.label,
            type=top.type,
            coverage_pct=top.coverage_pct,
            identity_pct=top.identity_pct,
            genes_present=sum(g.verdict == "present" for g in top.genes),
            genes_partial=sum(g.verdict == "partial" for g in top.genes),
            genes_absent=sum(g.verdict == "absent" for g in top.genes),
        )
    return ClusterReport(file="sample.fa", best=best, loci=tuple(loci))


class TestWeightedGenes:
    def test_all_present_scores_one(self) -> None:
        """Given every weighted gene present above the identity floor, When
        scored, Then the normalized score is exactly 1.0 with one positive
        normalized contribution per gene."""
        rule = WeightedGenesRule(
            model="weighted_genes",
            phenotype="K1",
            weights={"wzx": 1.0, "wzy": 3.0},
            identity_floor=95.0,
        )
        score, components = score_weighted_genes(rule, {"wzx": gene("wzx"), "wzy": gene("wzy")})
        assert score == 1.0
        assert components == (
            ScoreComponent(name="wzx", score=0.25),
            ScoreComponent(name="wzy", score=0.75),
        )

    def test_identity_floor_gates_presence(self) -> None:
        """Given a gene whose verdict is present but whose identity sits
        below the rule's floor, When scored, Then it does not count."""
        rule = WeightedGenesRule(
            model="weighted_genes", phenotype="K1", weights={"wzx": 1.0}, identity_floor=95.0
        )
        score, _ = score_weighted_genes(rule, {"wzx": gene("wzx", identity=94.99)})
        assert score == 0.0

    def test_partial_verdict_is_not_present(self) -> None:
        """Given a gene with the partial verdict, When scored, Then it does
        not count even at full identity."""
        rule = WeightedGenesRule(
            model="weighted_genes", phenotype="K1", weights={"wzx": 1.0}, identity_floor=90.0
        )
        score, _ = score_weighted_genes(rule, {"wzx": gene("wzx", verdict="partial")})
        assert score == 0.0

    def test_negative_gene_subtracts(self) -> None:
        """Given both positive genes and the negative gene present, When
        scored, Then the score is (1+1-1)/2 = 0.5."""
        rule = WeightedGenesRule(
            model="weighted_genes",
            phenotype="K1",
            weights={"a": 1.0, "b": 1.0},
            negative={"c": -1.0},
            identity_floor=90.0,
        )
        calls = {name: gene(name) for name in ("a", "b", "c")}
        score, components = score_weighted_genes(rule, calls)
        assert score == pytest.approx(0.5)
        assert components[-1] == ScoreComponent(name="c", score=-0.5)

    def test_absent_negative_gene_leaves_score_one(self) -> None:
        """Given the negative gene absent from the results, When scored,
        Then nothing is subtracted."""
        rule = WeightedGenesRule(
            model="weighted_genes",
            phenotype="K1",
            weights={"a": 1.0},
            negative={"c": -1.0},
            identity_floor=90.0,
        )
        score, _ = score_weighted_genes(rule, {"a": gene("a")})
        assert score == 1.0

    def test_require_any_unsatisfied_scores_zero(self) -> None:
        """Given a require_any set with no member present, When scored,
        Then the rule scores 0 regardless of the other weights."""
        rule = WeightedGenesRule(
            model="weighted_genes",
            phenotype="K1",
            weights={"a": 1.0},
            identity_floor=90.0,
            require_any=["wzx", "wzy"],
        )
        score, components = score_weighted_genes(rule, {"a": gene("a")})
        assert score == 0.0
        assert components == ()

    def test_require_any_satisfied_by_any_member(self) -> None:
        """Given one require_any member present, When scored, Then the rule
        fires normally."""
        rule = WeightedGenesRule(
            model="weighted_genes",
            phenotype="K1",
            weights={"a": 1.0},
            identity_floor=90.0,
            require_any=["wzx", "wzy"],
        )
        score, _ = score_weighted_genes(rule, {"a": gene("a"), "wzx": gene("wzx")})
        assert score == 1.0

    def test_denominator_zero_clamps_to_zero(self) -> None:
        """Given only negative weights with the negative gene present, When
        scored, Then the raw negative sum is clamped into [0, 1] as 0.0."""
        rule = WeightedGenesRule(
            model="weighted_genes",
            phenotype="K1",
            weights={},
            negative={"c": -1.0},
            identity_floor=90.0,
        )
        score, _ = score_weighted_genes(rule, {"c": gene("c")})
        assert score == 0.0


class TestClusterMatch:
    RULE = ClusterMatchRule(
        model="cluster_match",
        phenotype="K2",
        coverage=CoverageComponent(weight=0.5, floor=90.0),
        identity=IdentityComponent(weight=0.3, floor=95.0),
        key_genes=KeyGenesComponent(weight=0.2, genes=["wzx", "manC"]),
    )

    def test_perfect_locus_scores_weight_sum(self) -> None:
        """Given a locus at 100/100 with both key genes present, When
        scored, Then the total is the documented weight sum (1.0 here)."""
        genes = [gene("wzx"), gene("manC")]
        score, components = score_cluster_match(self.RULE, report([locus_call("locusA", genes)]))
        assert score == pytest.approx(1.0)
        assert [c.name for c in components] == ["coverage", "identity", "key_genes"]
        assert [c.score for c in components] == pytest.approx([0.5, 0.3, 0.2])

    def test_coverage_below_floor_scores_zero(self) -> None:
        """Given a locus below the coverage floor, When scored, Then the
        rule scores 0 for that locus."""
        genes = [gene("wzx"), gene("manC")]
        scored = report([locus_call("locusA", genes, coverage=89.9)])
        assert score_cluster_match(self.RULE, scored)[0] == 0.0

    def test_identity_below_floor_scores_zero(self) -> None:
        """Given a locus below the identity floor, When scored, Then the
        rule scores 0 even with perfect coverage."""
        genes = [gene("wzx"), gene("manC")]
        scored = report([locus_call("locusA", genes, identity=94.9)])
        assert score_cluster_match(self.RULE, scored)[0] == 0.0

    def test_key_gene_fraction_counts_present_only(self) -> None:
        """Given one of two key genes present, When scored, Then the
        key-genes component carries half its weight."""
        genes = [gene("wzx"), gene("manC", verdict="partial")]
        score, components = score_cluster_match(self.RULE, report([locus_call("locusA", genes)]))
        assert components[2].name == "key_genes"
        assert components[2].score == pytest.approx(0.1)
        assert score == pytest.approx(0.9)

    def test_partial_coverage_scales_component(self) -> None:
        """Given a locus at 95% coverage, When scored, Then the coverage
        component is weight * 0.95."""
        genes = [gene("wzx"), gene("manC")]
        _, components = score_cluster_match(
            self.RULE, report([locus_call("locusA", genes, coverage=95.0)])
        )
        assert components[0].name == "coverage"
        assert components[0].score == pytest.approx(0.475)

    def test_higher_scoring_locus_wins(self) -> None:
        """Given two covered loci, When scored, Then the higher total wins
        (here: the rank-2 locus whose key genes are present)."""
        weak = locus_call("locusA", [gene("gtrA")], coverage=100.0, rank=1)
        strong = locus_call("locusB", [gene("wzx"), gene("manC")], coverage=100.0, rank=2)
        score, components = score_cluster_match(self.RULE, report([weak, strong]))
        assert score == pytest.approx(1.0)
        assert [component.score for component in components] == pytest.approx([0.5, 0.3, 0.2])

    def test_no_loci_scores_zero(self) -> None:
        """Given a report with no covered loci, When scored, Then the rule
        scores 0 with no components."""
        assert score_cluster_match(self.RULE, report([])) == (0.0, ())

    def test_key_genes_absent_from_locus_count_as_absent(self) -> None:
        """Given a key gene that exists elsewhere but not in the scored
        locus, When scored, Then it counts as absent for this locus."""
        genes = [gene("wzx")]
        _, components = score_cluster_match(self.RULE, report([locus_call("locusA", genes)]))
        assert components[2].score == pytest.approx(0.1)


class TestLearnedLinear:
    def test_sigmoid_of_weighted_features_and_bias(self) -> None:
        """Given one present-gene feature with weight 2.0 and bias -1.0,
        When scored, Then the score is sigmoid(1)."""
        rule = LearnedLinearRule(
            model="learned_linear",
            phenotype="K3",
            features=["gene:wzx:present"],
            weights={"gene:wzx:present": 2.0},
            bias=-1.0,
            trained_on={"corpus": "fixture"},
        )
        score, components = score_learned_linear(
            rule, report([locus_call("l", [gene("wzx")])]), {"wzx": gene("wzx")}
        )
        assert score == pytest.approx(1.0 / (1.0 + math.exp(-1.0)))
        assert components == (
            ScoreComponent(name="bias", score=-1.0),
            ScoreComponent(name="gene:wzx:present", score=2.0),
        )

    def test_gene_cov_and_ident_fractions(self) -> None:
        """Given a gene at 80% coverage and 95% identity, When scored, Then
        cov/ident features carry 0.8/0.95 and present carries 1.0."""
        rule = LearnedLinearRule(
            model="learned_linear",
            phenotype="K3",
            features=["gene:wzx:cov", "gene:wzx:ident", "gene:wzx:present"],
            weights={"gene:wzx:cov": 1.0, "gene:wzx:ident": 1.0, "gene:wzx:present": -1.0},
            bias=0.0,
            trained_on={},
        )
        wzx = gene("wzx", coverage=80.0, identity=95.0)
        _, components = score_learned_linear(rule, report([locus_call("l", [wzx])]), {"wzx": wzx})
        assert [c.score for c in components[1:]] == pytest.approx([0.8, 0.95, -1.0])

    def test_cluster_features_read_locus_metrics(self) -> None:
        """Given cluster:<locus>:coverage/identity features, When scored,
        Then they read the named locus's fractions; an uncovered locus is 0."""
        rule = LearnedLinearRule(
            model="learned_linear",
            phenotype="K3",
            features=["cluster:locusA:coverage", "cluster:locusB:identity"],
            weights={"cluster:locusA:coverage": 1.0, "cluster:locusB:identity": 1.0},
            bias=0.0,
            trained_on={},
        )
        loci = [locus_call("locusA", [gene("wzx")], coverage=90.0, identity=99.0)]
        _, components = score_learned_linear(rule, report(loci), {})
        assert [c.score for c in components[1:]] == pytest.approx([0.9, 0.0])

    def test_missing_gene_feature_value_is_zero(self) -> None:
        """Given a gene feature for a gene with no call, When scored, Then
        the feature value is 0."""
        rule = LearnedLinearRule(
            model="learned_linear",
            phenotype="K3",
            features=["gene:absent_gene:present"],
            weights={"gene:absent_gene:present": 5.0},
            bias=0.0,
            trained_on={},
        )
        score, _ = score_learned_linear(
            rule, report([locus_call("l", [gene("wzx")])]), {"wzx": gene("wzx")}
        )
        assert score == pytest.approx(0.5)

    def test_feature_without_weight_is_malformed(self) -> None:
        """Given a feature missing from the weights map, When scored, Then
        a TYPING_MALFORMED DatabaseError names the feature."""
        rule = LearnedLinearRule(
            model="learned_linear",
            phenotype="K3",
            features=["gene:wzx:present"],
            weights={},
            bias=0.0,
            trained_on={},
        )
        with pytest.raises(DatabaseError) as raised:
            score_learned_linear(
                rule, report([locus_call("l", [gene("wzx")])]), {"wzx": gene("wzx")}
            )
        assert raised.value.code == "TYPING_MALFORMED"
        assert raised.value.context["feature"] == "gene:wzx:present"

    def test_malformed_feature_string_is_malformed(self) -> None:
        """Given a feature that is not gene:<id>:<metric> or
        cluster:<locus>:<metric>, When scored, Then TYPING_MALFORMED."""
        rule = LearnedLinearRule(
            model="learned_linear",
            phenotype="K3",
            features=["wzx"],
            weights={"wzx": 1.0},
            bias=0.0,
            trained_on={},
        )
        with pytest.raises(DatabaseError) as raised:
            score_learned_linear(rule, report([]), {})
        assert raised.value.code == "TYPING_MALFORMED"


def document(rules: list[object], cutoff: float = 0.9, margin: float = 0.05) -> TypingDocument:
    if not rules:  # the schema demands >= 1 rule; decision tests score separately
        rules = [
            WeightedGenesRule(
                model="weighted_genes", phenotype="unused", weights={}, identity_floor=0.0
            )
        ]
    return TypingDocument.model_validate(
        {
            "schema": "gapit.typing/1",
            "rules": rules,
            "cutoff": cutoff,
            "ambiguity_margin": margin,
            "fallback": "unknown",
        }
    )


class TestDecision:
    def test_clear_winner_is_called_with_high_confidence(self) -> None:
        """Given a top rule above cutoff separated from second by at least
        the margin, When decided, Then the phenotype is called with high
        confidence and the runner-up recorded."""
        doc = document([])
        phenotype, detail = decide_typing(
            [
                ScoredRule(phenotype="K1", score=0.95, components=()),
                ScoredRule(phenotype="K2", score=0.5, components=()),
            ],
            doc,
        )
        assert phenotype == "K1"
        assert detail.confidence == "high"
        assert detail.score == 0.95
        assert detail.runner_up is not None
        assert (detail.runner_up.phenotype, detail.runner_up.score) == ("K2", 0.5)
        assert detail.ambiguous == ()

    def test_close_scores_are_ambiguous(self) -> None:
        """Given a top rule above cutoff but separated by less than the
        margin, When decided, Then phenotype is None and the two tied
        phenotypes are listed."""
        phenotype, detail = decide_typing(
            [
                ScoredRule(phenotype="K1", score=0.95, components=()),
                ScoredRule(phenotype="K2", score=0.92, components=()),
            ],
            document([]),
        )
        assert phenotype is None
        assert detail.confidence == "ambiguous"
        assert [(entry.phenotype, entry.score) for entry in detail.ambiguous] == [
            ("K1", 0.95),
            ("K2", 0.92),
        ]

    def test_below_cutoff_falls_back_with_low_confidence(self) -> None:
        """Given every rule below cutoff, When decided, Then the fallback
        string is returned with low confidence."""
        phenotype, detail = decide_typing(
            [ScoredRule(phenotype="K1", score=0.5, components=())], document([])
        )
        assert phenotype == "unknown"
        assert detail.confidence == "low"

    def test_single_rule_above_cutoff_is_called(self) -> None:
        """Given exactly one rule above cutoff, When decided, Then it is
        called (no second rule exists to create ambiguity)."""
        phenotype, detail = decide_typing(
            [ScoredRule(phenotype="K1", score=0.95, components=())], document([])
        )
        assert phenotype == "K1"
        assert detail.confidence == "high"
        assert detail.runner_up is None

    def test_exact_tie_breaks_by_declared_order(self) -> None:
        """Given two rules with identical scores, When decided, Then the
        tie ranks by declared order (margin 0 lets an exact tie win)."""
        phenotype, detail = decide_typing(
            [
                ScoredRule(phenotype="K1", score=0.95, components=()),
                ScoredRule(phenotype="K2", score=0.95, components=()),
            ],
            document([], margin=0.0),
        )
        assert phenotype == "K1"
        assert detail.confidence == "high"
        assert detail.runner_up is not None and detail.runner_up.phenotype == "K2"


class TestEvaluateTyping:
    DOC = document(
        [
            WeightedGenesRule(
                model="weighted_genes",
                phenotype="K1",
                weights={"wzx": 1.0, "manC": 1.0},
                identity_floor=95.0,
                require_any=["wzx"],
            ),
            WeightedGenesRule(
                model="weighted_genes", phenotype="K2", weights={"orfX": 1.0}, identity_floor=95.0
            ),
        ]
    )

    def test_fills_phenotype_and_detail_on_best(self) -> None:
        """Given a report whose best locus carries wzx and manC present,
        When evaluated, Then best carries phenotype K1 with a detail whose
        components explain it."""
        loci = [locus_call("locusA", [gene("wzx"), gene("manC")])]
        typed = evaluate_typing(report(loci), self.DOC)
        assert typed.best is not None
        assert typed.best.phenotype == "K1"
        detail = typed.best.phenotype_detail
        assert detail is not None
        assert detail.confidence == "high"
        assert detail.score == 1.0
        assert detail.runner_up is not None and detail.runner_up.phenotype == "K2"
        assert [c.name for c in detail.components] == ["wzx", "manC"]

    def test_report_without_best_is_unchanged(self) -> None:
        """Given a report with no best call, When evaluated, Then the report
        is returned unchanged (nothing to annotate)."""
        low = report([locus_call("locusA", [gene("wzx")], coverage=50.0)])
        assert low.best is None
        assert evaluate_typing(low, self.DOC) == low

    def test_fallback_phenotype_when_below_cutoff(self) -> None:
        """Given no rule clearing the cutoff, When evaluated, Then the
        phenotype is the document's fallback string."""
        loci = [locus_call("locusA", [gene("wzx", verdict="absent")])]
        typed = evaluate_typing(report(loci), self.DOC)
        assert typed.best is not None
        assert typed.best.phenotype == "unknown"
        assert typed.best.phenotype_detail is not None
        assert typed.best.phenotype_detail.confidence == "low"


FEATURES = FeaturesDocument(
    loci=(
        LocusFeatures(
            id="locusA",
            label="A",
            type="KL101",
            genes=(GeneFeature(gene_id="wzx", start=1, end=100, strand="+"),),
        ),
        LocusFeatures(
            id="locusB",
            label="B",
            type="KL102",
            genes=(GeneFeature(gene_id="orfX", start=1, end=100, strand="+"),),
        ),
    )
)


class TestValidateReferences:
    def test_known_ids_pass(self) -> None:
        """Given rules referencing only ids the features carry, When
        validated, Then nothing raises."""
        validate_references(self._doc(), FEATURES)

    def test_unknown_gene_raises_typing_unknown_gene(self) -> None:
        """Given a rule referencing a gene no locus carries, When
        validated, Then TYPING_UNKNOWN_GENE names the gene and rule."""
        doc = document(
            [
                WeightedGenesRule(
                    model="weighted_genes",
                    phenotype="K1",
                    weights={"nope": 1.0},
                    identity_floor=90.0,
                )
            ]
        )
        with pytest.raises(DatabaseError) as raised:
            validate_references(doc, FEATURES)
        assert raised.value.code == "TYPING_UNKNOWN_GENE"
        assert raised.value.context["gene"] == "nope"
        assert raised.value.context["rule"] == "K1"

    def test_unknown_cluster_feature_locus_raises(self) -> None:
        """Given a learned feature naming a locus the db lacks, When
        validated, Then TYPING_UNKNOWN_GENE carries kind=locus."""
        doc = document(
            [
                LearnedLinearRule(
                    model="learned_linear",
                    phenotype="K3",
                    features=["cluster:locusZ:coverage"],
                    weights={"cluster:locusZ:coverage": 1.0},
                    bias=0.0,
                    trained_on={},
                )
            ]
        )
        with pytest.raises(DatabaseError) as raised:
            validate_references(doc, FEATURES)
        assert raised.value.code == "TYPING_UNKNOWN_GENE"
        assert raised.value.context["locus"] == "locusZ"

    def test_malformed_feature_raises_malformed(self) -> None:
        """Given a learned feature with a bad shape, When validated, Then
        TYPING_MALFORMED (before any id checking)."""
        doc = document(
            [
                LearnedLinearRule(
                    model="learned_linear",
                    phenotype="K3",
                    features=["bad"],
                    weights={"bad": 1.0},
                    bias=0.0,
                    trained_on={},
                )
            ]
        )
        with pytest.raises(DatabaseError) as raised:
            validate_references(doc, FEATURES)
        assert raised.value.code == "TYPING_MALFORMED"

    @staticmethod
    def _doc() -> TypingDocument:
        return document(
            [
                WeightedGenesRule(
                    model="weighted_genes",
                    phenotype="K1",
                    weights={"wzx": 1.0},
                    negative={"orfX": -0.5},
                    identity_floor=90.0,
                    require_any=["wzx"],
                ),
                ClusterMatchRule(
                    model="cluster_match",
                    phenotype="K2",
                    coverage=CoverageComponent(weight=0.5, floor=90.0),
                    identity=IdentityComponent(weight=0.5, floor=90.0),
                    key_genes=KeyGenesComponent(weight=0.0, genes=["wzx", "orfX"]),
                ),
                LearnedLinearRule(
                    model="learned_linear",
                    phenotype="K3",
                    features=["gene:wzx:present", "cluster:locusA:identity"],
                    weights={"gene:wzx:present": 1.0, "cluster:locusA:identity": 1.0},
                    bias=0.0,
                    trained_on={},
                ),
            ]
        )


class TestLoadTyping:
    def test_missing_file_returns_none(self, tmp_path: Path) -> None:
        """Given a cluster db directory without typing.json, When loaded,
        Then the result is None (untyped screening)."""
        database = Database(
            name="cps", path=tmp_path, sequences_path=tmp_path / "sequences", kind="cluster"
        )
        assert load_typing(database) is None

    def test_present_file_loads_document(self, tmp_path: Path) -> None:
        """Given typing.json in the db directory, When loaded, Then the
        validated gapit.typing/1 document comes back."""
        (tmp_path / "typing.json").write_text(
            TypingDocument.model_validate(
                {
                    "schema": "gapit.typing/1",
                    "rules": [
                        {
                            "model": "weighted_genes",
                            "phenotype": "K1",
                            "weights": {"wzx": 1.0},
                            "identity_floor": 95.0,
                        }
                    ],
                    "cutoff": 0.9,
                    "ambiguity_margin": 0.05,
                    "fallback": "unknown",
                }
            ).model_dump_json(),
            encoding="utf-8",
        )
        database = Database(
            name="cps", path=tmp_path, sequences_path=tmp_path / "sequences", kind="cluster"
        )
        loaded = load_typing(database)
        assert loaded is not None
        assert loaded.cutoff == 0.9


# ScoredRule must remain a plain value type (sorting + field access in the
# decision layer); pin its shape.
def test_scored_rule_is_a_value_tuple() -> None:
    """Given two ScoredRule values, When compared and indexed, Then they
    behave as named tuples (documented engine contract)."""
    rule = ScoredRule(phenotype="K1", score=0.5, components=(ScoreComponent(name="a", score=0.5),))
    assert rule.phenotype == "K1"
    assert rule.score == 0.5
    adapter = TypeAdapter(ScoredRule)
    assert adapter.validate_python(rule) == rule
