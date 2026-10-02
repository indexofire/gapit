"""Typing-engine stage-2 suites: the typing/2 rule primitives.

No I/O: the same direct ClusterReport/GeneCall values as the stage-1 engine
suites. Covers exact_set scoring + its documented tie exception,
coverage_floor gating, the scheme-level control-gene gate and unique-group
mixed-infection override, rule notes surfacing, and the pure compose
renderer. Schema-level validation of the same fields lives in
test_typing_models.py; the end-to-end fixture in test_gene_typing2.py.
"""

from typing import ClassVar, Literal

from gapit.cluster import GeneCall
from gapit.cluster_math import Verdict
from gapit.typing_decide import decide_scheme, decide_typing
from gapit.typing_engine import score_exact_set, score_weighted_genes
from gapit.typing_gene import compose_call
from gapit.typing_models import TypingScheme
from gapit.typing_results import PhenotypeDetail, PhenotypeScore, SchemeCall, ScoredRule
from gapit.typing_rules import ExactSetRule, WeightedGenesRule


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


def scheme(payload: dict[str, object]) -> TypingScheme:
    return TypingScheme.model_validate(payload)


def exact_rule(
    phenotype: str, requires: list[str], excludes: list[str] | None = None, **extra: object
) -> dict[str, object]:
    payload: dict[str, object] = {
        "model": "exact_set",
        "phenotype": phenotype,
        "requires": requires,
    }
    if excludes is not None:
        payload["excludes"] = excludes
    payload.update(extra)
    return payload


def weighted_rule(phenotype: str, weights: dict[str, float], **extra: object) -> dict[str, object]:
    payload: dict[str, object] = {
        "model": "weighted_genes",
        "phenotype": phenotype,
        "weights": weights,
        "identity_floor": 90.0,
    }
    payload.update(extra)
    return payload


class TestExactSetScoring:
    def test_satisfied_rule_scores_one(self) -> None:
        """Given every required gene present and every excluded gene absent,
        When scored, Then the score is exactly 1.0 with one 1.0 component
        per gene (required then excluded)."""
        rule = ExactSetRule.model_validate(
            exact_rule("4b-4d-4e", ["orf2110", "orf2870"], ["lmo1118"])
        )
        score, components = score_exact_set(
            rule, {"orf2110": gene("orf2110"), "orf2870": gene("orf2870")}
        )
        assert score == 1.0
        assert [(c.name, c.score) for c in components] == [
            ("orf2110", 1.0),
            ("orf2870", 1.0),
            ("lmo1118", 1.0),
        ]

    def test_missing_required_gene_scores_zero(self) -> None:
        """Given one required gene with no call at all, When scored, Then
        the rule scores 0.0 and that gene's component is 0.0."""
        rule = ExactSetRule.model_validate(exact_rule("1/2a-3a", ["lmo0737", "prs"]))
        score, components = score_exact_set(rule, {"lmo0737": gene("lmo0737")})
        assert score == 0.0
        assert [(c.name, c.score) for c in components] == [("lmo0737", 1.0), ("prs", 0.0)]

    def test_present_excluded_gene_scores_zero(self) -> None:
        """Given an excluded gene present above both floors, When scored,
        Then the rule scores 0.0 with that gene's component 0.0."""
        rule = ExactSetRule.model_validate(exact_rule("1/2c-3c", ["lmo1118"], ["lmo0737"]))
        score, components = score_exact_set(
            rule, {"lmo1118": gene("lmo1118"), "lmo0737": gene("lmo0737")}
        )
        assert score == 0.0
        assert (components[-1].name, components[-1].score) == ("lmo0737", 0.0)

    def test_default_floors_are_ninety_ninety(self) -> None:
        """Given a required gene at 90/90 exactly versus 89.9 identity and
        89.9 coverage, When scored without explicit floors, Then only the
        90/90 gene counts (the documented defaults)."""
        rule = ExactSetRule.model_validate(exact_rule("X", ["g"]))
        assert score_exact_set(rule, {"g": gene("g", identity=90.0, coverage=90.0)})[0] == 1.0
        assert score_exact_set(rule, {"g": gene("g", identity=89.9)})[0] == 0.0
        assert score_exact_set(rule, {"g": gene("g", coverage=89.9)})[0] == 0.0

    def test_partial_verdict_never_counts(self) -> None:
        """Given a required gene with the partial verdict, When scored,
        Then the rule scores 0 even at full identity/coverage."""
        rule = ExactSetRule.model_validate(exact_rule("X", ["g"]))
        assert score_exact_set(rule, {"g": gene("g", verdict="partial")})[0] == 0.0


class TestExactSetRequiresAny:
    def test_any_of_satisfied_by_one_present_gene(self) -> None:
        """Given a requires_any set with exactly one member present (no
        plain requires), When scored, Then the rule scores 1.0 with each
        any-of gene's presence as its component."""
        rule = ExactSetRule.model_validate(
            exact_rule("ETEC", [], requires_any=["lt", "sth", "stp"])
        )
        score, components = score_exact_set(rule, {"sth": gene("sth")})
        assert score == 1.0
        assert [(c.name, c.score) for c in components] == [
            ("lt", 0.0),
            ("sth", 1.0),
            ("stp", 0.0),
        ]

    def test_any_of_unsatisfied_scores_zero(self) -> None:
        """Given a requires_any set with no member present — genes missing
        entirely or below a floor — When scored, Then the rule scores
        0.0 (the any-of gate, like weighted_genes require_any)."""
        rule = ExactSetRule.model_validate(exact_rule("EAEC", [], requires_any=["aggR", "pic"]))
        assert score_exact_set(rule, {})[0] == 0.0
        below_floor = score_exact_set(rule, {"aggR": gene("aggR", identity=89.9)})[0]
        assert below_floor == 0.0

    def test_combined_all_of_any_of_and_none_of(self) -> None:
        """Given a rule combining requires, requires_any, and excludes,
        When scored, Then 1.0 iff all three constraints hold, with the
        components listing requires, then requires_any, then excludes."""
        rule = ExactSetRule.model_validate(
            exact_rule(
                "EHEC",
                ["escV"],
                ["bfpB"],
                requires_any=["stx1a", "stx2a"],
            )
        )
        score, components = score_exact_set(
            rule,
            {"escV": gene("escV"), "stx2a": gene("stx2a")},
        )
        assert score == 1.0
        assert [(c.name, c.score) for c in components] == [
            ("escV", 1.0),
            ("stx1a", 0.0),
            ("stx2a", 1.0),
            ("bfpB", 1.0),
        ]
        no_any_of = score_exact_set(rule, {"escV": gene("escV")})[0]
        broken_all_of = score_exact_set(rule, {"stx2a": gene("stx2a")})[0]
        present_exclude = score_exact_set(
            rule, {"escV": gene("escV"), "stx2a": gene("stx2a"), "bfpB": gene("bfpB")}
        )[0]
        assert (no_any_of, broken_all_of, present_exclude) == (0.0, 0.0, 0.0)

    def test_absent_requires_any_keeps_requires_only_rules_working(self) -> None:
        """Given a rule that omits requires_any entirely (the pre-existing
        shape), When scored, Then only requires/excludes decide (backward
        compatibility: the field is optional and inert when empty)."""
        rule = ExactSetRule.model_validate(exact_rule("EIEC", ["invE"]))
        assert score_exact_set(rule, {"invE": gene("invE")})[0] == 1.0
        assert score_exact_set(rule, {})[0] == 0.0


class TestExactSetTieSemantics:
    def test_satisfied_exact_set_wins_tie_by_declaration_order(self) -> None:
        """Given a satisfied exact_set at 1.0 declared before a
        weighted_genes rule also at 1.0, When decided with a nonzero margin,
        Then the call is the exact_set's phenotype with high confidence
        (the documented exception: ordering decides, not the margin)."""
        payload = scheme(
            {
                "name": "doumith",
                "rules": [
                    exact_rule("1/2a-3a", ["lmo0737"]),
                    weighted_rule("other", {"lmo0737": 1.0}),
                ],
                "cutoff": 0.9,
                "ambiguity_margin": 0.05,
                "fallback": "NT",
            }
        )
        scored = [
            ScoredRule("1/2a-3a", 1.0, (), exact=True),
            ScoredRule("other", 1.0, ()),
        ]
        phenotype, detail = decide_typing(scored, payload)
        assert phenotype == "1/2a-3a"
        assert detail.confidence == "high"

    def test_exact_set_declared_second_loses_the_same_tie(self) -> None:
        """Given the same 1.0/1.0 tie with the weighted rule declared
        first, When decided, Then the weighted rule wins — the exception is
        symmetric ordering, not exact_set priority."""
        payload = scheme(
            {
                "name": "doumith",
                "rules": [
                    weighted_rule("other", {"lmo0737": 1.0}),
                    exact_rule("1/2a-3a", ["lmo0737"]),
                ],
                "cutoff": 0.9,
                "ambiguity_margin": 0.05,
                "fallback": "NT",
            }
        )
        scored = [
            ScoredRule("other", 1.0, ()),
            ScoredRule("1/2a-3a", 1.0, (), exact=True),
        ]
        phenotype, _ = decide_typing(scored, payload)
        assert phenotype == "other"

    def test_unsatisfied_exact_set_does_not_claim_the_exception(self) -> None:
        """Given an exact_set at 0.9 versus a weighted rule at 0.88 with
        margin 0.05, When decided, Then the call stays ambiguous-by-margin
        (the exception applies only to a satisfied 1.0 tie)."""
        payload = scheme(
            {
                "name": "doumith",
                "rules": [
                    exact_rule("a", ["g"]),
                    weighted_rule("b", {"g": 1.0}),
                ],
                "cutoff": 0.5,
                "ambiguity_margin": 0.05,
                "fallback": "NT",
            }
        )
        phenotype, detail = decide_typing(
            [ScoredRule("a", 0.9, (), exact=True), ScoredRule("b", 0.88, ())], payload
        )
        assert phenotype is None
        assert detail.confidence == "ambiguous"

    def test_non_exact_ties_stay_ambiguous(self) -> None:
        """Given two weighted rules tied at 1.0 with margin 0.05, When
        decided, Then the call is ambiguous (stage-1 semantics untouched)."""
        payload = scheme(
            {
                "name": "k",
                "rules": [
                    weighted_rule("K13", {"wzy_a": 1.0}),
                    weighted_rule("K64", {"wzy_b": 1.0}),
                ],
                "cutoff": 0.9,
                "ambiguity_margin": 0.05,
                "fallback": "K?",
            }
        )
        phenotype, detail = decide_typing(
            [ScoredRule("K13", 1.0, ()), ScoredRule("K64", 1.0, ())], payload
        )
        assert phenotype is None
        assert detail.confidence == "ambiguous"


class TestCoverageFloor:
    def test_weighted_gene_below_coverage_floor_does_not_count(self) -> None:
        """Given a weighted_genes rule with coverage_floor 95 and a gene at
        94.9% coverage, When scored, Then the gene does not count; at 95 it
        does."""
        rule = WeightedGenesRule.model_validate(
            weighted_rule("K1", {"wzx": 1.0}, coverage_floor=95.0)
        )
        assert score_weighted_genes(rule, {"wzx": gene("wzx", coverage=94.9)})[0] == 0.0
        assert score_weighted_genes(rule, {"wzx": gene("wzx", coverage=95.0)})[0] == 1.0

    def test_weighted_genes_have_no_coverage_gate_by_default(self) -> None:
        """Given a present-verdict gene at only 50% coverage and a rule
        without coverage_floor, When scored, Then the gene still counts
        (no gate beyond the engine's verdict semantics)."""
        rule = WeightedGenesRule.model_validate(weighted_rule("K1", {"wzx": 1.0}))
        assert score_weighted_genes(rule, {"wzx": gene("wzx", coverage=50.0)})[0] == 1.0


class TestControlGene:
    SCHEME: ClassVar[dict[str, object]] = {
        "name": "doumith",
        "rules": [exact_rule("1/2a-3a", ["lmo0737"])],
        "cutoff": 0.9,
        "ambiguity_margin": 0.05,
        "fallback": "NT",
        "control_gene": "prs",
    }

    def test_absent_control_gene_zeroes_the_scheme(self) -> None:
        """Given the marker gene present but the control gene absent, When
        decided, Then the scheme outputs its fallback with low confidence,
        score 0.0, and the 'control gene absent' note."""
        phenotype, detail = decide_scheme(
            scheme(self.SCHEME),
            [ScoredRule("1/2a-3a", 1.0, (), exact=True)],
            {"lmo0737": gene("lmo0737")},
        )
        assert phenotype == "NT"
        assert detail == PhenotypeDetail(
            score=0.0, confidence="low", components=(), notes=("control gene absent",)
        )

    def test_partial_control_gene_is_absent(self) -> None:
        """Given the control gene with a partial verdict, When decided,
        Then the gate still fires (verdict semantics, not hit presence)."""
        phenotype, detail = decide_scheme(
            scheme(self.SCHEME),
            [ScoredRule("1/2a-3a", 1.0, (), exact=True)],
            {"prs": gene("prs", verdict="partial")},
        )
        assert phenotype == "NT"
        assert detail.confidence == "low"

    def test_present_control_gene_lets_the_scheme_decide(self) -> None:
        """Given the control gene present, When decided, Then the normal
        ranked decision runs and the marker rule is called."""
        phenotype, detail = decide_scheme(
            scheme(self.SCHEME),
            [ScoredRule("1/2a-3a", 1.0, (), exact=True)],
            {"prs": gene("prs"), "lmo0737": gene("lmo0737")},
        )
        assert phenotype == "1/2a-3a"
        assert detail.confidence == "high"


class TestUniqueGroup:
    SCHEME: ClassVar[dict[str, object]] = {
        "name": "o_group",
        "rules": [
            exact_rule("O1", ["wzx_o1"]),
            exact_rule("O2", ["wzx_o2"]),
        ],
        "cutoff": 0.9,
        "ambiguity_margin": 0.05,
        "fallback": "O?",
        "unique_group": {"wzx": ["wzx_o1", "wzx_o2"]},
        "mixed_phenotype": "mixed",
    }

    def test_two_present_group_genes_call_mixed(self) -> None:
        """Given both genes of the wzx unique group present, When decided,
        Then the phenotype is the mixed phenotype with the pair listed in
        ambiguous and a mixed-infection note (low confidence)."""
        phenotype, detail = decide_scheme(
            scheme(self.SCHEME),
            [ScoredRule("O1", 1.0, (), exact=True), ScoredRule("O2", 1.0, (), exact=True)],
            {"wzx_o1": gene("wzx_o1"), "wzx_o2": gene("wzx_o2")},
        )
        assert phenotype == "mixed"
        assert detail.confidence == "low"
        assert detail.score == 1.0
        assert [(entry.phenotype, entry.score) for entry in detail.ambiguous] == [
            ("wzx_o1", 1.0),
            ("wzx_o2", 1.0),
        ]
        assert detail.notes == ("mixed infection in unique group 'wzx'",)

    def test_single_present_group_gene_calls_normally(self) -> None:
        """Given only one member of the group present, When decided, Then
        the normal decision runs and that gene's rule is called."""
        phenotype, _ = decide_scheme(
            scheme(self.SCHEME),
            [ScoredRule("O1", 1.0, (), exact=True), ScoredRule("O2", 0.0, ())],
            {"wzx_o1": gene("wzx_o1")},
        )
        assert phenotype == "O1"

    def test_absent_group_genes_call_the_fallback(self) -> None:
        """Given no group member present, When decided, Then the group does
        not fire and the scheme falls back below cutoff."""
        phenotype, detail = decide_scheme(
            scheme(self.SCHEME), [ScoredRule("O1", 0.0, ()), ScoredRule("O2", 0.0, ())], {}
        )
        assert phenotype == "O?"
        assert detail.confidence == "low"


class TestRuleNotes:
    def test_winning_rule_notes_surface_verbatim(self) -> None:
        """Given a winning rule carrying notes, When decided, Then the
        detail surfaces them verbatim as a tuple."""
        payload = scheme(
            {
                "name": "doumith",
                "rules": [
                    exact_rule(
                        "4b-4d-4e",
                        ["orf2110"],
                        notes=["4b/4d/4e share one PCR pattern", "pINV-ipaB flag"],
                    )
                ],
                "cutoff": 0.9,
                "ambiguity_margin": 0.05,
                "fallback": "NT",
            }
        )
        _, detail = decide_typing(
            [ScoredRule("4b-4d-4e", 1.0, (), ("4b/4d/4e share one PCR pattern", "pINV-ipaB flag"))],
            payload,
        )
        assert detail.notes == ("4b/4d/4e share one PCR pattern", "pINV-ipaB flag")

    def test_notes_empty_when_rule_has_none(self) -> None:
        """Given a winning rule without notes, When decided, Then the
        detail carries no notes (the default)."""
        payload = scheme(
            {
                "name": "doumith",
                "rules": [exact_rule("X", ["g"])],
                "cutoff": 0.9,
                "ambiguity_margin": 0.05,
                "fallback": "NT",
            }
        )
        _, detail = decide_typing([ScoredRule("X", 1.0, ())], payload)
        assert detail.notes == ()


class TestCompose:
    def test_compose_renders_called_phenotypes(self) -> None:
        """Given a {o}:{k} template and both ingredients called with high
        confidence, When composed, Then the phenotype is the rendered
        string with the minimum ingredient score and one component per
        ingredient."""
        call = compose_call(
            "{o_group}:{k_group}",
            {
                "o_group": _scheme_call("O1", 1.0, "high"),
                "k_group": _scheme_call("K13", 0.75, "high"),
            },
        )
        assert call.phenotype == "O1:K13"
        assert call.confidence == "high"
        assert call.score == 0.75
        assert [(c.name, c.score) for c in call.components] == [("o_group", 1.0), ("k_group", 0.75)]

    def test_ingredient_fallback_flows_through_as_low_confidence(self) -> None:
        """Given one ingredient on its fallback string (low confidence),
        When composed, Then the fallback string renders into the value and
        the composed confidence is low."""
        call = compose_call(
            "{o_group}:{k_group}",
            {
                "o_group": _scheme_call("O?", 0.0, "low"),
                "k_group": _scheme_call("K13", 1.0, "high"),
            },
        )
        assert call.phenotype == "O?:K13"
        assert call.confidence == "low"

    def test_ambiguous_ingredient_yields_null_with_variants(self) -> None:
        """Given an ambiguous ingredient (phenotype null with variants),
        When composed, Then the composed phenotype is null, confidence
        ambiguous, and the ambiguous list carries the ingredient's
        variants."""
        call = compose_call(
            "{o_group}:{k_group}",
            {
                "o_group": _scheme_call("O1", 1.0, "high"),
                "k_group": _scheme_call(None, 1.0, "ambiguous", variants=("K13", "K64")),
            },
        )
        assert call.phenotype is None
        assert call.confidence == "ambiguous"
        assert [entry.phenotype for entry in call.ambiguous] == ["K13", "K64"]
        assert call.notes == ("ingredient scheme 'k_group' is ambiguous",)

    def test_mixed_ingredient_string_flows_through(self) -> None:
        """Given an ingredient whose mixed-infection override produced a
        non-null phenotype, When composed, Then that string renders like
        any other called value."""
        call = compose_call(
            "{o_group}:{k_group}",
            {
                "o_group": _scheme_call("mixed", 1.0, "low"),
                "k_group": _scheme_call("K13", 1.0, "high"),
            },
        )
        assert call.phenotype == "mixed:K13"
        assert call.confidence == "low"


def _scheme_call(
    phenotype: str | None,
    score: float,
    confidence: Literal["high", "ambiguous", "low"],
    variants: tuple[str, ...] = (),
) -> SchemeCall:
    return SchemeCall(
        phenotype=phenotype,
        score=score,
        confidence=confidence,
        ambiguous=tuple(PhenotypeScore(phenotype=name, score=1.0) for name in variants),
    )
