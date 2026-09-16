"""Unit tests for BLAST outfmt-6 row parsing and the blastn version gate."""

import pytest

from gapit.blast import BLAST_FIELDS, BlastRow, ensure_blast, parse_blast_row
from gapit.errors import GapitError

# Captured verbatim from the real pipeline (any2fasta | blastn, BLAST+ 2.17).
CANONICAL = (
    "contig1\t1\t79\t79\ttinyamr~~~tetA~~~NC_000913.3:100-900~~~TETRACYCLINE\t"
    "1\t79\t79\tplus\t8.72e-40\t79\t100.000\t0\t0\t"
    "tinyamr~~~tetA~~~NC_000913.3:100-900~~~TETRACYCLINE tetracycline efflux pump TetA"
)


def test_blast_fields_exact_order() -> None:
    """The 15 outfmt-6 fields in SPEC.md §3 order."""
    assert BLAST_FIELDS == [
        "qseqid",
        "qstart",
        "qend",
        "qlen",
        "sseqid",
        "sstart",
        "send",
        "slen",
        "sstrand",
        "evalue",
        "length",
        "pident",
        "gaps",
        "gapopen",
        "stitle",
    ]


def test_parse_canonical_row() -> None:
    """Given a real 15-field line, When parsed, Then every field lands typed."""
    parsed = parse_blast_row(CANONICAL)
    assert isinstance(parsed, BlastRow)
    assert parsed.qseqid == "contig1"
    assert (parsed.qstart, parsed.qend, parsed.qlen) == (1, 79, 79)
    assert parsed.sseqid == "tinyamr~~~tetA~~~NC_000913.3:100-900~~~TETRACYCLINE"
    assert (parsed.sstart, parsed.send, parsed.slen) == (1, 79, 79)
    assert parsed.sstrand == "plus"
    assert parsed.evalue == pytest.approx(8.72e-40)
    assert parsed.length == 79
    assert parsed.pident == 100.0
    assert (parsed.gaps, parsed.gapopen) == (0, 0)
    assert parsed.stitle.endswith("tetracycline efflux pump TetA")


def test_parse_zero_evalue() -> None:
    """Given evalue '0.0' (ultra-significant hits), When parsed, Then it is a float."""
    line = CANONICAL.replace("\t8.72e-40\t", "\t0.0\t")
    assert parse_blast_row(line).evalue == 0.0


def test_row_with_14_fields_fails() -> None:
    """Given a line missing the stitle column, When parsed, Then BLAST_PARSE_FAILED
    with the upstream wording."""
    fields = CANONICAL.split("\t")
    with pytest.raises(GapitError) as excinfo:
        parse_blast_row("\t".join(fields[:14]))
    assert excinfo.value.code == "BLAST_PARSE_FAILED"
    assert str(excinfo.value) == "can not find sequence data"


def test_row_with_16_fields_fails() -> None:
    """Given a line with an extra column, When parsed, Then BLAST_PARSE_FAILED."""
    with pytest.raises(GapitError) as excinfo:
        parse_blast_row(CANONICAL + "\textra")
    assert excinfo.value.code == "BLAST_PARSE_FAILED"


def test_ensure_blast_passes_version_gate() -> None:
    """Given the env's blastn (>= 2.7), When ensured, Then no exception."""
    assert ensure_blast() is None
