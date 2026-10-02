"""Tests for the gapit.typing schema (gapit.typing_models) — versions 1 and 2.

Schema + validation only — no evaluation logic. All three rule model types
must validate; malformed documents must surface as typed DatabaseErrors with
context. gapit.typing/2 adds the named multi-scheme container; gapit.typing/1
(flat rules/cutoff/margin/fallback) must degrade to one anonymous ``default``
scheme so existing typed cluster databases keep their behavior.
"""

import json
from pathlib import Path

import pytest

from gapit.errors import DatabaseError, InputError
from gapit.gbfeatures import FeaturesDocument, GeneFeature, LocusFeatures
from gapit.typing_gene import validate_gene_typing
from gapit.typing_models import (
    V1_SCHEME_NAME,
    TypingDocument,
    TypingDocumentV1,
    read_typing_document,
    single_scheme,
    template_placeholders,
    typing_schema_of,
    validate_references,
)
from gapit.typing_rules import ClusterMatchRule, LearnedLinearRule, WeightedGenesRule

DATA = Path(__file__).parent / "data" / "cluster"

WEIGHTED_GENES_RULE = {
    "model": "weighted_genes",
    "phenotype": "K6",
    "weights": {"wzx": 1.0, "manC": 2.0},
    "negative": {"orfA": -1.0},
    "identity_floor": 95.0,
    "require_any": ["wzx", "manC"],
}


def v2_document(schemes: list[dict[str, object]]) -> dict[str, object]:
    return {"schema": "gapit.typing/2", "schemes": schemes}


def scheme_document(name: str = "serotype") -> dict[str, object]:
    return {
        "name": name,
        "rules": [WEIGHTED_GENES_RULE],
        "cutoff": 0.9,
        "ambiguity_margin": 0.05,
        "fallback": "unknown",
    }


def write_document(tmp_path: Path, payload: dict[str, object], filename: str = "doc.json") -> Path:
    path = tmp_path / filename
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


class TestV1Degradation:
    def test_valid_document_roundtrips_all_three_model_types(self) -> None:
        """Given a typing/1 document with one rule of each model type, When
        read, Then it loads as ONE default scheme whose rules validate with
        their model kinds, per-model fields carry through, and the
        scheme-level cutoff/margin/fallback survive a JSON round-trip."""
        document = read_typing_document(DATA / "typing_valid.json")

        assert document.schema_name == "gapit.typing/2"
        assert len(document.schemes) == 1
        scheme = document.schemes[0]
        assert scheme.name == V1_SCHEME_NAME
        assert [rule.model for rule in scheme.rules] == [
            "weighted_genes",
            "cluster_match",
            "learned_linear",
        ]
        weighted, cluster, learned = scheme.rules
        assert isinstance(weighted, WeightedGenesRule)
        assert isinstance(cluster, ClusterMatchRule)
        assert isinstance(learned, LearnedLinearRule)
        assert weighted.phenotype == "K6"
        assert weighted.weights == {"wzx": 1.0, "manC": 2.0}
        assert weighted.negative == {"orfA": -1.0}
        assert weighted.identity_floor == 95.0
        assert weighted.require_any == ["wzx", "manC"]
        assert cluster.phenotype == "K7"
        assert (cluster.coverage.weight, cluster.coverage.floor) == (0.6, 90.0)
        assert (cluster.identity.weight, cluster.identity.floor) == (0.4, 95.0)
        assert (cluster.key_genes.weight, cluster.key_genes.genes) == (1.0, ["wzx", "manC"])
        assert learned.phenotype == "K8"
        assert learned.features == ["gene:wzx:cov", "gene:manC:present", "cluster:KL9001:coverage"]
        assert learned.bias == -0.5
        assert learned.trained_on == {"corpus": "vp-serotypes-2026", "version": "1"}
        assert scheme.cutoff == 0.9
        assert scheme.ambiguity_margin == 0.05
        assert scheme.fallback == "unknown"
        assert json.loads(document.model_dump_json(by_alias=True))["schemes"][0]["rules"][0][
            "model"
        ] == ("weighted_genes")

    def test_v1_document_equals_its_v2_mirror(self) -> None:
        """Given the same spec as a typing/1 file and as a typing/2 file
        with one scheme, When both are read, Then the parsed documents are
        equal (name aside, the v1 scheme name IS the mirror's name)."""
        v1 = read_typing_document(DATA / "typing_valid.json")
        v1_direct = TypingDocumentV1.model_validate_json(
            (DATA / "typing_valid.json").read_text(encoding="utf-8")
        )
        mirror = v1_direct.as_document()
        assert v1 == mirror
        assert single_scheme(v1).rules == v1_direct.rules

    def test_unknown_model_type_is_typed_database_error(self) -> None:
        """Given a rule whose model is not one of the three kinds, When read,
        Then a DatabaseError TYPING_MALFORMED names the file."""
        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(DATA / "typing_invalid.json")

        assert excinfo.value.code == "TYPING_MALFORMED"
        assert excinfo.value.context["file"].endswith("typing_invalid.json")
        assert "neural_net" in str(excinfo.value)

    def test_negative_floor_is_rejected(self, tmp_path: Path) -> None:
        """Given a weighted_genes rule with a negative identity floor, When read,
        Then a typed DatabaseError rejects the document."""
        bad = write_document(
            tmp_path,
            {
                "schema": "gapit.typing/1",
                "rules": [
                    {
                        "model": "weighted_genes",
                        "phenotype": "K6",
                        "weights": {"wzx": 1.0},
                        "identity_floor": -5.0,
                    }
                ],
                "cutoff": 0.9,
                "ambiguity_margin": 0.05,
                "fallback": "unknown",
            },
        )

        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(bad)

        assert excinfo.value.code == "TYPING_MALFORMED"

    def test_missing_cutoff_is_rejected(self, tmp_path: Path) -> None:
        """Given a structurally valid rules list but no document-level cutoff,
        When read, Then a typed DatabaseError rejects the document."""
        bad = write_document(
            tmp_path,
            {
                "schema": "gapit.typing/1",
                "rules": [
                    {
                        "model": "weighted_genes",
                        "phenotype": "K6",
                        "weights": {"wzx": 1.0},
                        "identity_floor": 90.0,
                    }
                ],
                "ambiguity_margin": 0.05,
                "fallback": "unknown",
            },
        )

        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(bad)

        assert excinfo.value.code == "TYPING_MALFORMED"

    def test_missing_file_is_input_not_found(self, tmp_path: Path) -> None:
        """Given a --typing path that does not exist, When read, Then the typed
        INPUT_NOT_FOUND error propagates (the --tsv precedent: a missing
        user-supplied file is an input error, exit 5, not a db error)."""
        with pytest.raises(InputError) as excinfo:
            read_typing_document(tmp_path / "absent.json")

        assert excinfo.value.code == "INPUT_NOT_FOUND"
        assert excinfo.value.context["file"].endswith("absent.json")


class TestV2Documents:
    def test_multi_scheme_document_parses(self, tmp_path: Path) -> None:
        """Given a typing/2 document with two named schemes, When read, Then
        both schemes parse with their own rules and thresholds, in declared
        order."""
        path = write_document(
            tmp_path,
            v2_document(
                [
                    scheme_document("pathotype"),
                    {
                        "name": "toxin",
                        "rules": [
                            {
                                "model": "weighted_genes",
                                "phenotype": "Stx2",
                                "weights": {"stx2": 1.0},
                                "identity_floor": 90.0,
                            }
                        ],
                        "cutoff": 0.95,
                        "ambiguity_margin": 0.1,
                        "fallback": "none",
                    },
                ]
            ),
        )

        document = read_typing_document(path)

        assert document.schema_name == "gapit.typing/2"
        assert [scheme.name for scheme in document.schemes] == ["pathotype", "toxin"]
        assert document.schemes[1].cutoff == 0.95
        assert document.schemes[1].fallback == "none"
        assert document.schemes[0].rules[0].phenotype == "K6"

    def test_duplicate_scheme_names_are_rejected(self, tmp_path: Path) -> None:
        """Given two schemes sharing a name, When read, Then TYPING_MALFORMED
        names the duplicate."""
        path = write_document(
            tmp_path, v2_document([scheme_document("same"), scheme_document("same")])
        )

        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(path)

        assert excinfo.value.code == "TYPING_MALFORMED"
        assert "duplicate scheme names" in str(excinfo.value)
        assert "same" in str(excinfo.value)

    @pytest.mark.parametrize("bad_name", ["Pathotype", "patho type", "patho-1", "1pathotype!", ""])
    def test_scheme_name_charset_is_enforced(self, tmp_path: Path, bad_name: str) -> None:
        """Given a scheme name outside [a-z0-9_], When read, Then
        TYPING_MALFORMED rejects the document."""
        path = write_document(tmp_path, v2_document([scheme_document(bad_name)]))

        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(path)

        assert excinfo.value.code == "TYPING_MALFORMED"

    def test_missing_scheme_name_is_rejected(self, tmp_path: Path) -> None:
        """Given a v2 scheme without a name (report/1 keys phenotypes by
        name), When read, Then TYPING_MALFORMED."""
        scheme = scheme_document()
        del scheme["name"]
        path = write_document(tmp_path, v2_document([scheme]))

        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(path)

        assert excinfo.value.code == "TYPING_MALFORMED"

    def test_empty_schemes_list_is_rejected(self, tmp_path: Path) -> None:
        """Given a v2 document with no schemes, When read, Then
        TYPING_MALFORMED."""
        path = write_document(tmp_path, v2_document([]))

        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(path)

        assert excinfo.value.code == "TYPING_MALFORMED"

    def test_not_json_is_malformed(self, tmp_path: Path) -> None:
        """Given a typing file that is not JSON at all, When read, Then
        TYPING_MALFORMED (not a raw json error)."""
        path = tmp_path / "bad.json"
        path.write_text("{not json", encoding="utf-8")

        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(path)

        assert excinfo.value.code == "TYPING_MALFORMED"


class TestSchemaTag:
    def test_v1_and_v2_tags_are_reported(self, tmp_path: Path) -> None:
        """Given installed v1 and v2 documents, When tagged, Then the literal
        schema string comes back (the manifest's typing_schema value)."""
        v1_path = write_document(
            tmp_path,
            {
                "schema": "gapit.typing/1",
                "rules": [WEIGHTED_GENES_RULE],
                "cutoff": 0.9,
                "ambiguity_margin": 0.05,
                "fallback": "unknown",
            },
            filename="v1.json",
        )
        v2_path = write_document(tmp_path, v2_document([scheme_document()]), filename="v2.json")

        assert typing_schema_of(v1_path) == "gapit.typing/1"
        assert typing_schema_of(v2_path) == "gapit.typing/2"

    def test_unknown_tag_is_malformed(self, tmp_path: Path) -> None:
        """Given a document with an unknown schema tag, When tagged, Then
        TYPING_MALFORMED."""
        path = write_document(tmp_path, {"schema": "gapit.typing/9"})

        with pytest.raises(DatabaseError) as excinfo:
            typing_schema_of(path)

        assert excinfo.value.code == "TYPING_MALFORMED"


class TestSingleScheme:
    def test_multi_scheme_document_is_rejected_for_cluster_use(self, tmp_path: Path) -> None:
        """Given a two-scheme document and the cluster path's one-phenotype
        slot, When single_scheme is called, Then TYPING_MALFORMED names both
        scheme names."""
        path = write_document(
            tmp_path, v2_document([scheme_document("pathotype"), scheme_document("toxin")])
        )
        document = read_typing_document(path)

        with pytest.raises(DatabaseError) as excinfo:
            single_scheme(document)

        assert excinfo.value.code == "TYPING_MALFORMED"
        assert "pathotype" in str(excinfo.value) and "toxin" in str(excinfo.value)

    def test_single_scheme_document_returns_its_scheme(self, tmp_path: Path) -> None:
        """Given a one-scheme document (any version), When single_scheme is
        called, Then that scheme comes back."""
        path = write_document(tmp_path, v2_document([scheme_document("serotype")]))

        assert single_scheme(read_typing_document(path)).name == "serotype"
        assert single_scheme(read_typing_document(DATA / "typing_valid.json")).name == (
            V1_SCHEME_NAME
        )


EXACT_SET_RULE = {
    "model": "exact_set",
    "phenotype": "4b-4d-4e",
    "requires": ["orf2110", "orf2870"],
    "notes": ["4b/4d/4e share one PCR pattern"],
}

FEATURES = FeaturesDocument(
    loci=(
        LocusFeatures(
            id="locusA",
            label="A",
            type="KL101",
            genes=(GeneFeature(gene_id="prs", start=1, end=100, strand="+"),),
        ),
    )
)


class TestStage2RuleFields:
    def test_exact_set_rule_roundtrips(self, tmp_path: Path) -> None:
        """Given an exact_set rule with floors and notes, When read, Then
        every field survives (defaults 90/90 when floors are omitted)."""
        path = write_document(
            tmp_path,
            v2_document([scheme_document("doumith") | {"rules": [EXACT_SET_RULE]}]),
        )

        scheme = read_typing_document(path).schemes[0]

        (rule,) = scheme.rules
        assert rule.model == "exact_set"
        assert rule.requires == ["orf2110", "orf2870"]
        assert rule.excludes == []
        assert rule.identity_floor == 90.0
        assert rule.coverage_floor == 90.0
        assert rule.notes == ["4b/4d/4e share one PCR pattern"]

    def test_exact_set_overlap_is_malformed(self, tmp_path: Path) -> None:
        """Given an exact_set rule requiring AND excluding the same gene,
        When read, Then TYPING_MALFORMED."""
        rule = EXACT_SET_RULE | {"excludes": ["orf2110"]}
        path = write_document(
            tmp_path, v2_document([scheme_document("doumith") | {"rules": [rule]}])
        )

        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(path)

        assert excinfo.value.code == "TYPING_MALFORMED"
        assert "both required and excluded" in str(excinfo.value)

    def test_exact_set_without_requires_is_malformed(self, tmp_path: Path) -> None:
        """Given an exact_set rule with an empty requires list, When read,
        Then TYPING_MALFORMED (a vacuous always-1.0 rule)."""
        rule = EXACT_SET_RULE | {"requires": []}
        path = write_document(
            tmp_path, v2_document([scheme_document("doumith") | {"rules": [rule]}])
        )

        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(path)

        assert excinfo.value.code == "TYPING_MALFORMED"

    def test_exact_set_requires_any_roundtrips(self, tmp_path: Path) -> None:
        """Given an exact_set rule carrying only a requires_any set (the
        DEC any-of shape), When read, Then it validates with requires
        defaulting to empty and every field surviving."""
        rule = EXACT_SET_RULE | {"requires": [], "requires_any": ["stx1a", "stx2a"]}
        path = write_document(tmp_path, v2_document([scheme_document("dec") | {"rules": [rule]}]))

        parsed = read_typing_document(path).schemes[0].rules[0]

        assert parsed.model == "exact_set"
        assert parsed.requires == []
        assert parsed.requires_any == ["stx1a", "stx2a"]

    @pytest.mark.parametrize(
        ("overlap_rule", "fragment"),
        [
            ({"requires_any": ["orf2110"]}, "required and in requires_any"),
            (
                {"requires": [], "requires_any": ["orf2110"], "excludes": ["orf2110"]},
                "in requires_any and excluded",
            ),
        ],
    )
    def test_exact_set_requires_any_overlap_is_malformed(
        self, tmp_path: Path, overlap_rule: dict[str, object], fragment: str
    ) -> None:
        """Given an exact_set rule whose requires_any overlaps requires or
        excludes, When read, Then TYPING_MALFORMED naming the conflict."""
        rule = EXACT_SET_RULE | overlap_rule
        path = write_document(
            tmp_path, v2_document([scheme_document("doumith") | {"rules": [rule]}])
        )

        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(path)

        assert excinfo.value.code == "TYPING_MALFORMED"
        assert fragment in str(excinfo.value)

    @pytest.mark.parametrize("floor", [-0.1, 100.1])
    def test_out_of_range_coverage_floor_is_malformed(self, tmp_path: Path, floor: float) -> None:
        """Given a coverage floor outside [0, 100] on either rule kind,
        When read, Then TYPING_MALFORMED."""
        exact = EXACT_SET_RULE | {"coverage_floor": floor}
        weighted = WEIGHTED_GENES_RULE | {"coverage_floor": floor}
        for rules in ([exact], [weighted]):
            path = write_document(
                tmp_path, v2_document([scheme_document("doumith") | {"rules": rules}])
            )
            with pytest.raises(DatabaseError) as excinfo:
                read_typing_document(path)
            assert excinfo.value.code == "TYPING_MALFORMED"


class TestStage2SchemeFields:
    def test_control_unique_mixed_fields_roundtrip(self, tmp_path: Path) -> None:
        """Given a scheme with control_gene, unique_group, and
        mixed_phenotype, When read, Then all three survive together with
        the exact_set rules."""
        scheme = scheme_document("o_group") | {
            "rules": [
                {"model": "exact_set", "phenotype": "O1", "requires": ["wzx_o1"]},
                {"model": "exact_set", "phenotype": "O2", "requires": ["wzx_o2"]},
            ],
            "control_gene": "prs",
            "unique_group": {"wzx": ["wzx_o1", "wzx_o2"]},
            "mixed_phenotype": "mixed",
        }
        path = write_document(tmp_path, v2_document([scheme]))

        parsed = read_typing_document(path).schemes[0]

        assert parsed.control_gene == "prs"
        assert parsed.unique_group == {"wzx": ["wzx_o1", "wzx_o2"]}
        assert parsed.mixed_phenotype == "mixed"

    def test_scheme_without_rules_or_compose_is_malformed(self, tmp_path: Path) -> None:
        """Given a scheme declaring neither rules nor compose, When read,
        Then TYPING_MALFORMED."""
        scheme = scheme_document("empty")
        scheme["rules"] = []
        path = write_document(tmp_path, v2_document([scheme]))

        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(path)

        assert excinfo.value.code == "TYPING_MALFORMED"

    @pytest.mark.parametrize(
        ("drop", "fragment"),
        [("mixed", "mixed_phenotype"), ("group", "unique_group")],
    )
    def test_unique_group_and_mixed_need_each_other(
        self, tmp_path: Path, drop: str, fragment: str
    ) -> None:
        """Given unique_group without mixed_phenotype (or the reverse),
        When read, Then TYPING_MALFORMED (the pair must be set together)."""
        scheme: dict[str, object] = {
            "name": "o_group",
            "rules": [{"model": "exact_set", "phenotype": "O1", "requires": ["wzx_o1"]}],
            "cutoff": 0.9,
            "ambiguity_margin": 0.05,
            "fallback": "O?",
            "unique_group": {"wzx": ["wzx_o1", "wzx_o2"]},
            "mixed_phenotype": "mixed",
        }
        del scheme[fragment]
        path = write_document(tmp_path, v2_document([scheme]))

        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(path)

        assert excinfo.value.code == "TYPING_MALFORMED"
        assert "together" in str(excinfo.value)

    def test_compose_scheme_with_rules_is_malformed(self, tmp_path: Path) -> None:
        """Given a compose scheme that also declares rules, When read,
        Then TYPING_MALFORMED."""
        scheme = scheme_document("serotype") | {"compose": "{serotype}:{serotype}"}
        path = write_document(tmp_path, v2_document([scheme]))

        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(path)

        assert excinfo.value.code == "TYPING_MALFORMED"
        assert "compose scheme declares no rules" in str(excinfo.value)


class TestComposeValidation:
    def ingredients(self) -> list[dict[str, object]]:
        return [
            {
                "name": "o_group",
                "rules": [{"model": "exact_set", "phenotype": "O1", "requires": ["wzx_o1"]}],
                "cutoff": 0.9,
                "ambiguity_margin": 0.05,
                "fallback": "O?",
            },
            {
                "name": "k_group",
                "rules": [{"model": "exact_set", "phenotype": "K13", "requires": ["wzy_k13"]}],
                "cutoff": 0.9,
                "ambiguity_margin": 0.05,
                "fallback": "K?",
            },
        ]

    def compose_scheme(self, template: str) -> dict[str, object]:
        return {
            "name": "serotype",
            "compose": template,
            "cutoff": 0.9,
            "ambiguity_margin": 0.05,
            "fallback": "O?:K?",
        }

    def test_valid_compose_document_parses(self, tmp_path: Path) -> None:
        """Given a compose scheme over two sibling schemes, When read, Then
        the document parses with the template intact."""
        path = write_document(
            tmp_path, v2_document([*self.ingredients(), self.compose_scheme("{o_group}:{k_group}")])
        )

        document = read_typing_document(path)

        assert document.schemes[2].compose == "{o_group}:{k_group}"

    @pytest.mark.parametrize(
        "template",
        ["{o_group}:{missing}", "{serotype}", "{o_group", "{}", "{o_group:>5}"],
    )
    def test_bad_templates_are_malformed(self, tmp_path: Path, template: str) -> None:
        """Given a compose template naming an unknown scheme, itself, or
        using non-plain brace syntax, When read, Then TYPING_MALFORMED."""
        path = write_document(
            tmp_path, v2_document([*self.ingredients(), self.compose_scheme(template)])
        )

        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(path)

        assert excinfo.value.code == "TYPING_MALFORMED"

    def test_placeholder_naming_a_compose_scheme_is_malformed(self, tmp_path: Path) -> None:
        """Given a compose scheme whose placeholder names another compose
        scheme, When read, Then TYPING_MALFORMED (compose composes called
        schemes only — no ordering chains)."""
        schemes = [
            *self.ingredients(),
            self.compose_scheme("{o_group}:{k_group}"),
            self.compose_scheme("{o_group}:{serotype}") | {"name": "banner"},
        ]
        path = write_document(tmp_path, v2_document(schemes))

        with pytest.raises(DatabaseError) as excinfo:
            read_typing_document(path)

        assert excinfo.value.code == "TYPING_MALFORMED"

    def test_template_placeholders_dedupe_in_order(self) -> None:
        """Given a template repeating a placeholder, When extracted, Then
        the names come back deduplicated in first-occurrence order."""
        assert template_placeholders("{b}:{a}:{b}") == ("b", "a")


class TestStage2References:
    def document_with(self, scheme: dict[str, object], tmp_path: Path) -> TypingDocument:
        return read_typing_document(write_document(tmp_path, v2_document([scheme])))

    def test_unknown_exact_set_gene_raises(self, tmp_path: Path) -> None:
        """Given an exact_set rule referencing a gene the features lack,
        When validated, Then TYPING_UNKNOWN_GENE names the gene and rule."""
        scheme = scheme_document("doumith") | {"rules": [EXACT_SET_RULE]}
        with pytest.raises(DatabaseError) as raised:
            validate_references(self.document_with(scheme, tmp_path), FEATURES)
        assert raised.value.code == "TYPING_UNKNOWN_GENE"
        assert raised.value.context["gene"] == "orf2110"
        assert raised.value.context["rule"] == "4b-4d-4e"

    def test_unknown_requires_any_gene_raises(self, tmp_path: Path) -> None:
        """Given an exact_set rule whose requires_any references a gene the
        features lack, When validated, Then TYPING_UNKNOWN_GENE (the any-of
        set is checked like every other gene reference)."""
        rule = {
            "model": "exact_set",
            "phenotype": "X",
            "requires": ["prs"],
            "requires_any": ["aggR"],
        }
        scheme = scheme_document("dec") | {"rules": [rule]}
        with pytest.raises(DatabaseError) as raised:
            validate_references(self.document_with(scheme, tmp_path), FEATURES)
        assert raised.value.code == "TYPING_UNKNOWN_GENE"
        assert raised.value.context["gene"] == "aggR"

    def test_unknown_control_gene_raises(self, tmp_path: Path) -> None:
        """Given a scheme whose control gene is not in the features, When
        validated, Then TYPING_UNKNOWN_GENE names the scheme."""
        scheme = scheme_document("doumith") | {
            "rules": [{"model": "exact_set", "phenotype": "X", "requires": ["prs"]}],
            "control_gene": "ipaH",
        }
        with pytest.raises(DatabaseError) as raised:
            validate_references(self.document_with(scheme, tmp_path), FEATURES)
        assert raised.value.code == "TYPING_UNKNOWN_GENE"
        assert raised.value.context["scheme"] == "doumith"
        assert raised.value.context["gene"] == "ipaH"

    def test_unknown_unique_group_member_raises(self, tmp_path: Path) -> None:
        """Given a unique group with a member the features lack, When
        validated, Then TYPING_UNKNOWN_GENE names the scheme."""
        scheme = scheme_document("o_group") | {
            "rules": [{"model": "exact_set", "phenotype": "O1", "requires": ["prs"]}],
            "unique_group": {"wzx": ["prs", "wzx_o1"]},
            "mixed_phenotype": "mixed",
        }
        with pytest.raises(DatabaseError) as raised:
            validate_references(self.document_with(scheme, tmp_path), FEATURES)
        assert raised.value.code == "TYPING_UNKNOWN_GENE"
        assert raised.value.context["gene"] == "wzx_o1"

    def test_gene_build_checks_scheme_genes(self, tmp_path: Path) -> None:
        """Given a gene-path document whose unique group references a gene
        the FASTA lacks, When the build validates, Then
        TYPING_UNKNOWN_GENE names the scheme."""
        scheme = scheme_document("o_group") | {
            "rules": [{"model": "exact_set", "phenotype": "O1", "requires": ["prs"]}],
            "unique_group": {"wzx": ["prs", "wzx_o1"]},
            "mixed_phenotype": "mixed",
        }
        with pytest.raises(DatabaseError) as raised:
            validate_gene_typing(self.document_with(scheme, tmp_path), frozenset({"prs"}))
        assert raised.value.code == "TYPING_UNKNOWN_GENE"
        assert raised.value.context["scheme"] == "o_group"
