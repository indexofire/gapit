"""gapit.serotype/1 rule evaluation tests."""

from __future__ import annotations

import pytest

from gspec.cluster import GeneCall
from gspec.serotype_engine import evaluate_serotype
from gspec.serotype_models import SerotypeDocument, SerotypeRule


def _call(
    gene: str,
    *,
    coverage: float = 100.0,
    identity: float = 100.0,
    alignment_start: float | None = None,
) -> GeneCall:
    return GeneCall(
        gene_id=gene,
        start=1,
        end=100,
        strand="+",
        coverage_pct=coverage,
        identity_pct=identity,
        verdict="present",
        alignment_start_pct=alignment_start,
    )


def _document() -> SerotypeDocument:
    return SerotypeDocument(
        name="test",
        control_gene="control",
        rules=(
            SerotypeRule(phenotype="type A", requires=("a",)),
            SerotypeRule(phenotype="type B", requires=("b",)),
            SerotypeRule(
                phenotype="exempt type",
                requires=("x",),
                control_exempt=True,
            ),
        ),
    )


def test_serotype_rule_call_is_high_confidence() -> None:
    call = evaluate_serotype({"control": _call("control"), "a": _call("a")}, _document())
    assert (call.phenotype, call.confidence) == ("type A", "high")


def test_multiple_matching_rules_are_ambiguous() -> None:
    call = evaluate_serotype(
        {"control": _call("control"), "a": _call("a"), "b": _call("b")},
        _document(),
    )
    assert call.phenotype is None
    assert call.confidence == "ambiguous"
    assert {entry.phenotype for entry in call.ambiguous} == {"type A", "type B"}


def test_missing_control_falls_back_without_control_exempt() -> None:
    call = evaluate_serotype({"a": _call("a")}, _document())
    assert (call.phenotype, call.confidence) == ("Not Shigella/EIEC", "low")


def test_control_exempt_rule_fires_without_control() -> None:
    call = evaluate_serotype({"x": _call("x")}, _document())
    assert (call.phenotype, call.confidence) == ("exempt type", "high")


def test_rule_floors_filter_weak_markers() -> None:
    document = SerotypeDocument(
        name="test",
        rules=(
            SerotypeRule(
                phenotype="strict",
                requires=("a",),
                coverage_floor=90.0,
                identity_floor=95.0,
            ),
        ),
    )
    call = evaluate_serotype({"a": _call("a", coverage=80.0)}, document)
    assert call.phenotype == "Shigella spp. (unresolved)"


def _sb6_sb10_document() -> SerotypeDocument:
    """The shipped SB6/SB10 rule shape: a 5'-intact wbaM is type 10, a
    5'-truncated (IS629) wbaM is type 6, and missing wbaM stays ambiguous."""
    return SerotypeDocument(
        name="test",
        rules=(
            SerotypeRule(
                phenotype="S. boydii 10",
                requires=("Sb6_wzx", "wbaM"),
                max_alignment_start_pct={"wbaM": 20.0},
            ),
            SerotypeRule(
                phenotype="S. boydii 6",
                requires=("Sb6_wzx", "wbaM"),
                max_coverage={"wbaM": 95.0},
                min_alignment_start_pct={"wbaM": 20.0},
            ),
            SerotypeRule(
                phenotype="S. boydii 6/10",
                requires=("Sb6_wzx",),
                excludes=("wbaM",),
            ),
        ),
    )


def test_wbam_alignment_start_separates_partial_assemblies() -> None:
    """Given wbaM alignments, When the 5' end is covered, Then the call is
    SB10 even with a 3'-truncated assembly; a 5'-truncated alignment (the
    IS629 insertion) is SB6."""
    document = _sb6_sb10_document()
    calls = {
        "Sb6_wzx": _call("Sb6_wzx"),
        "wbaM": _call("wbaM", alignment_start=0.0),
    }
    intact = evaluate_serotype(calls, document)
    three_prime_partial = evaluate_serotype(
        {"Sb6_wzx": _call("Sb6_wzx"), "wbaM": _call("wbaM", coverage=86.83, alignment_start=0.0)},
        document,
    )
    five_prime_partial = evaluate_serotype(
        {
            "Sb6_wzx": _call("Sb6_wzx"),
            "wbaM": _call("wbaM", coverage=73.67, alignment_start=26.33),
        },
        document,
    )
    assert intact.phenotype == "S. boydii 10"
    assert three_prime_partial.phenotype == "S. boydii 10"
    assert five_prime_partial.phenotype == "S. boydii 6"


def test_wbam_start_gate_requires_alignment_evidence() -> None:
    """Given reads tables (no alignment start), When wbaM is present, Then
    the SB6/SB10 rules do not match: breadth cannot locate the IS629
    junction, and the resolver stays unresolved instead of guessing."""
    document = _sb6_sb10_document()
    for coverage in (73.67, 100.0):
        call = evaluate_serotype(
            {"Sb6_wzx": _call("Sb6_wzx"), "wbaM": _call("wbaM", coverage=coverage)},
            document,
        )
        assert call.phenotype == "Shigella spp. (unresolved)"


def test_missing_wbam_keeps_sb6_sb10_ambiguous() -> None:
    """Given Sb6_wzx without wbaM, When the resolver runs, Then it keeps
    the explicit 6/10 ambiguity rather than guessing."""
    call = evaluate_serotype({"Sb6_wzx": _call("Sb6_wzx")}, _sb6_sb10_document())
    assert call.phenotype == "S. boydii 6/10"


def test_alignment_start_gates_validate_their_markers() -> None:
    """Given a gate on a non-positive marker or a non-percentage value,
    When the rule is constructed, Then validation rejects it."""
    with pytest.raises(ValueError, match="not a positive marker"):
        SerotypeRule(
            phenotype="bad",
            requires=("a",),
            max_alignment_start_pct={"other": 10.0},
        )
    with pytest.raises(ValueError, match="must be a percentage"):
        SerotypeRule(
            phenotype="bad",
            requires=("a",),
            min_alignment_start_pct={"a": 150.0},
        )
