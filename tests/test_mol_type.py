"""Unit tests for the abricate mol_type heuristic (SPEC.md §2).

Rule (replicate exactly): delete [AGTC] case-insensitively; if the remaining
length is STRICTLY greater than 50% of the original, the db is protein.
"""

from gaita.db import mol_type


def test_pure_acgt_is_nucl() -> None:
    assert mol_type("ATGGTCAATTCCGCTGGCGCTGGCGATGGTTCGGCAACCGTGCGCTGGCA") == "nucl"


def test_lowercase_acgt_counts_as_nucleotide() -> None:
    assert mol_type("agtcagtcag") == "nucl"


def test_clearly_protein_letters_are_prot() -> None:
    assert mol_type("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAG") == "prot"


def test_exactly_half_non_acgt_is_nucl() -> None:
    """6 letters, 2 non-ACGT: 2 > 3 is false -> nucl (strictly-greater rule)."""
    assert mol_type("AGTCNN") == "nucl"


def test_even_length_exactly_half_is_nucl() -> None:
    """8 letters, 4 non-ACGT: 4 > 4 is false -> nucl."""
    assert mol_type("AGTCNNNN") == "nucl"


def test_just_over_half_is_prot() -> None:
    """9 letters, 5 non-ACGT: 5 > 4.5 -> prot."""
    assert mol_type("AGTCNNNNN") == "prot"


def test_empty_string_is_nucl() -> None:
    assert mol_type("") == "nucl"
