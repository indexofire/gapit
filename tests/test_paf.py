"""Unit tests for PAF row parsing (SPEC.md §10)."""

import pytest

from gapit.errors import GapitError
from gapit.reads import parse_paf_row

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
