"""Unit tests for PAF row parsing (SPEC.md §10)."""

import pytest

from gapit.errors import GapitError
from gapit.paf import PafRecord, parse_paf_row

CANONICAL = (
    "r1\t150\t5\t60\t+\tdb~~~geneA~~~ACC~~~RES\t100\t10\t55\t50\t50\t60\ttp:A:P\tcm:i:12\ts1:i:50"
)


def test_parse_canonical_row_with_primary_tag() -> None:
    """Given a 12-field PAF line plus tags, When parsed, Then every field lands
    typed and tp:A:P marks the record primary."""
    record = parse_paf_row(CANONICAL)
    assert record.qname == "r1"
    assert (record.qlen, record.qstart, record.qend) == (150, 5, 60)
    assert record.strand == "+"
    assert record.tname == "db~~~geneA~~~ACC~~~RES"
    assert (record.tlen, record.tstart, record.tend) == (100, 10, 55)
    assert (record.nmatch, record.alen, record.mapq) == (50, 50, 60)
    assert record.is_primary is True


def test_parse_secondary_tag() -> None:
    """Given tp:A:S, When parsed, Then the record is not primary."""
    record = parse_paf_row(CANONICAL.replace("tp:A:P", "tp:A:S"))
    assert record.is_primary is False


def test_missing_tp_tag_counts_as_primary() -> None:
    """Given a row without any tp tag, When parsed, Then it is primary."""
    bare = CANONICAL.split("\t")[:12]
    assert parse_paf_row("\t".join(bare)).is_primary is True


def test_row_with_too_few_fields_fails() -> None:
    """Given an 11-field line, When parsed, Then PAF_PARSE_FAILED."""
    fields = CANONICAL.split("\t")[:11]
    with pytest.raises(GapitError) as excinfo:
        parse_paf_row("\t".join(fields))
    assert excinfo.value.code == "PAF_PARSE_FAILED"


def test_parse_minus_strand() -> None:
    record = parse_paf_row(CANONICAL.replace("\t+\t", "\t-\t"))
    assert record.strand == "-"


def test_parse_equals_fully_validated_record() -> None:
    """Given the canonical row (and its tagged/untagged variants), When parsed
    via parse_paf_row (model_construct path), Then the result equals a fully
    validated PafRecord built from the same typed values — the manual
    int()/strand checks make pydantic re-validation redundant, not absent."""
    expected = PafRecord(
        qname="r1",
        qlen=150,
        qstart=5,
        qend=60,
        strand="+",
        tname="db~~~geneA~~~ACC~~~RES",
        tlen=100,
        tstart=10,
        tend=55,
        nmatch=50,
        alen=50,
        mapq=60,
        is_primary=True,
    )
    assert parse_paf_row(CANONICAL) == expected
    secondary = parse_paf_row(CANONICAL.replace("tp:A:P", "tp:A:S"))
    assert secondary == expected.model_copy(update={"is_primary": False})
    bare = "\t".join(CANONICAL.split("\t")[:12])
    assert parse_paf_row(bare) == expected
