"""Tests for the gapit.typing/1 declarative scoring schema (gapit.typing_models).

Schema + validation only — no evaluation logic in this stage. All three rule
model types must validate; malformed documents must surface as typed
DatabaseErrors with context.
"""

import json
from pathlib import Path

import pytest

from gapit.errors import DatabaseError, InputError
from gapit.typing_models import (
    ClusterMatchRule,
    LearnedLinearRule,
    WeightedGenesRule,
    read_typing_document,
)

DATA = Path(__file__).parent / "data" / "cluster"


def test_valid_document_roundtrips_all_three_model_types() -> None:
    """Given a typing document with one rule of each model type, When read,
    Then every rule validates with its model kind, per-model fields carry
    through, and the document-level cutoff/margin/fallback survive a JSON
    round-trip."""
    document = read_typing_document(DATA / "typing_valid.json")

    assert document.schema_name == "gapit.typing/1"
    assert [rule.model for rule in document.rules] == [
        "weighted_genes",
        "cluster_match",
        "learned_linear",
    ]
    weighted, cluster, learned = document.rules
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
    assert document.cutoff == 0.9
    assert document.ambiguity_margin == 0.05
    assert document.fallback == "unknown"
    assert json.loads(document.model_dump_json(by_alias=True))["rules"][0]["model"] == (
        "weighted_genes"
    )


def test_unknown_model_type_is_typed_database_error() -> None:
    """Given a rule whose model is not one of the three kinds, When read,
    Then a DatabaseError TYPING_MALFORMED names the file."""
    with pytest.raises(DatabaseError) as excinfo:
        read_typing_document(DATA / "typing_invalid.json")

    assert excinfo.value.code == "TYPING_MALFORMED"
    assert excinfo.value.context["file"].endswith("typing_invalid.json")
    assert "neural_net" in str(excinfo.value)


def test_negative_floor_is_rejected(tmp_path: Path) -> None:
    """Given a weighted_genes rule with a negative identity floor, When read,
    Then a typed DatabaseError rejects the document."""
    bad = tmp_path / "neg.json"
    bad.write_text(
        json.dumps(
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
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(DatabaseError) as excinfo:
        read_typing_document(bad)

    assert excinfo.value.code == "TYPING_MALFORMED"


def test_missing_cutoff_is_rejected(tmp_path: Path) -> None:
    """Given a structurally valid rules list but no document-level cutoff,
    When read, Then a typed DatabaseError rejects the document."""
    bad = tmp_path / "nocut.json"
    bad.write_text(
        json.dumps(
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
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(DatabaseError) as excinfo:
        read_typing_document(bad)

    assert excinfo.value.code == "TYPING_MALFORMED"


def test_missing_file_is_input_not_found(tmp_path: Path) -> None:
    """Given a --typing path that does not exist, When read, Then the typed
    INPUT_NOT_FOUND error propagates (the --tsv precedent: a missing
    user-supplied file is an input error, exit 5, not a db error)."""
    with pytest.raises(InputError) as excinfo:
        read_typing_document(tmp_path / "absent.json")

    assert excinfo.value.code == "INPUT_NOT_FOUND"
    assert excinfo.value.context["file"].endswith("absent.json")
