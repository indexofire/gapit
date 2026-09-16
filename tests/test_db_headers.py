"""Unit tests for ~~~ database header parsing (SPEC.md §4 step 5)."""

from gapit.db import parse_db_header


def test_canonical_four_fields() -> None:
    """Given a full db~~~gene~~~accession~~~resistance header, When parsed,
    Then every field is populated."""
    header = parse_db_header("ncbi~~~tetA~~~NC_000913.3~~~TETRACYCLINE", default_db="ncbi")
    assert header.database == "ncbi"
    assert header.gene == "tetA"
    assert header.accession == "NC_000913.3"
    assert header.resistance == "TETRACYCLINE"


def test_no_separator_falls_back_to_gene_only() -> None:
    """Given a seqid with no ~~~ at all, When parsed, Then gene is the whole id,
    database is the default, and trailing fields are empty (1-field case)."""
    header = parse_db_header("tetA(1)", default_db="resfinder")
    assert header.gene == "tetA(1)"
    assert header.database == "resfinder"
    assert header.accession == ""
    assert header.resistance == ""


def test_two_fields_pad_trailing() -> None:
    """Given db~~~gene, When parsed, Then accession and resistance are empty."""
    header = parse_db_header("vfdb~~~toxA", default_db="ncbi")
    assert header.database == "vfdb"
    assert header.gene == "toxA"
    assert header.accession == ""
    assert header.resistance == ""


def test_three_fields_pad_resistance() -> None:
    """Given db~~~gene~~~accession, When parsed, Then resistance is empty."""
    header = parse_db_header("card~~~blaTEM-1~~~J01749.1", default_db="ncbi")
    assert header.database == "card"
    assert header.gene == "blaTEM-1"
    assert header.accession == "J01749.1"
    assert header.resistance == ""


def test_empty_database_field_uses_default() -> None:
    """Given a leading ~~~ (database field absent), When parsed, Then the default
    database is substituted."""
    header = parse_db_header("~~~gene~~~acc~~~RES", default_db="megares")
    assert header.database == "megares"
    assert header.gene == "gene"
    assert header.accession == "acc"
    assert header.resistance == "RES"


def test_resistance_classes_preserved_verbatim() -> None:
    """Given a multi-class ;-joined resistance field, When parsed, Then it is
    preserved verbatim."""
    header = parse_db_header(
        "argannot~~~gene~~~acc~~~BETA-LACTAM;AMINOGLYCOSIDE;MACROLIDE", default_db="ncbi"
    )
    assert header.resistance == "BETA-LACTAM;AMINOGLYCOSIDE;MACROLIDE"
