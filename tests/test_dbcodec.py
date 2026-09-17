"""Unit tests for the gapit/v1 tagged header codec (dbcodec, Phase 7 Wave A1)."""

import pytest

from gapit.db import parse_db_header
from gapit.dbcodec import decode_seqid, encode_seqid, is_gapit_header
from gapit.errors import DatabaseError


def test_encode_produces_tagged_wire_format() -> None:
    """Given db/gene/accession/function values containing escapable chars, When
    encoded, Then the wire format is gapit|db=..|gene=..|acc=..|func=.. with
    lowercase-hex percent escapes and fixed key order."""
    seqid = encode_seqid("ncbi", "te|t%A=1", "NC_1;2", ["TET", "AMP"])
    assert seqid == "gapit|db=ncbi|gene=te%7Ct%25A%3D1|acc=NC_1%3B2|func=TET%3BAMP"


def test_roundtrip_plain_probe() -> None:
    """Given a plain alphanumeric header (mixed case, dashes, dots), When
    encoded then decoded, Then every field round-trips exactly."""
    header = decode_seqid(
        encode_seqid("card", "blaTEM-1a", "J01749.1", ["BETA-LACTAM"]), default_db="ncbi"
    )
    assert header.database == "card"
    assert header.gene == "blaTEM-1a"
    assert header.accession == "J01749.1"
    assert header.function == "BETA-LACTAM"


def test_roundtrip_probe_escape_zoo() -> None:
    """Given values containing every probe-verified escapable char (| = % ; space),
    When encoded then decoded, Then all values survive byte-exact."""
    header = decode_seqid(
        encode_seqid("re|sfinder%", "ge|ne%=x;y z", "ac|c%=1", ["a|b%=c;d e"]),
        default_db="ncbi",
    )
    assert header.database == "re|sfinder%"
    assert header.gene == "ge|ne%=x;y z"
    assert header.accession == "ac|c%=1"
    assert header.function == "a|b%=c;d e"


def test_roundtrip_empty_function() -> None:
    """Given an empty function list, When encoded then decoded, Then function
    is the empty string."""
    header = decode_seqid(encode_seqid("ncbi", "tetA", "NC_1", []), default_db="ncbi")
    assert header.function == ""


def test_roundtrip_multi_function_display_form() -> None:
    """Given a multi-class function list, When encoded then decoded, Then
    function is the abricate display string 'a;b;c'."""
    header = decode_seqid(
        encode_seqid("argannot", "gene", "acc", ["BETA-LACTAM", "AMINOGLYCOSIDE", "MACROLIDE"]),
        default_db="ncbi",
    )
    assert header.function == "BETA-LACTAM;AMINOGLYCOSIDE;MACROLIDE"


def test_roundtrip_727_char_header() -> None:
    """Given a header whose total seqid length exceeds the 727-char probe size,
    When encoded then decoded, Then the gene round-trips without truncation."""
    gene = "x" * 700
    seqid = encode_seqid("ncbi", gene, "acc", ["TET"])
    assert len(seqid) >= 727
    assert decode_seqid(seqid, default_db="ncbi").gene == gene


def test_encode_is_total_for_weird_input() -> None:
    """Given values with newlines, tabs, CR, unicode and all-empty fields, When
    encoded, Then it never raises and weird values round-trip exactly."""
    weird = "β-lacta\tmase\n100%|raw=|=;\r\n🌊"
    seqid = encode_seqid("", weird, weird, [weird, weird])
    assert is_gapit_header(seqid)
    assert decode_seqid(seqid, default_db="vfdb").gene == weird
    assert decode_seqid(seqid, default_db="vfdb").function == f"{weird};{weird}"
    assert encode_seqid("", "", "", []) == "gapit|db=|gene=|acc=|func="


def test_decode_skips_unknown_keys() -> None:
    """Given a header with an extra unknown key (forward compat), When decoded,
    Then the known fields parse and the unknown key is skipped."""
    header = decode_seqid("gapit|db=ncbi|gene=tetA|future=x|acc=NC_1", default_db="ncbi")
    assert header.database == "ncbi"
    assert header.gene == "tetA"
    assert header.accession == "NC_1"


def test_decode_accepts_uppercase_hex_escapes() -> None:
    """Given a value escaped with uppercase hex, When decoded, Then it decodes
    (encode emits lowercase, decode reads both cases)."""
    assert decode_seqid("gapit|db=ncbi|gene=a%7Cb|acc=NC_1", default_db="ncbi").gene == "a|b"


def test_decode_db_fallback() -> None:
    """Given a new-format header whose db key is present but empty, When decoded,
    Then database falls back to default_db (symmetric with legacy)."""
    empty = decode_seqid("gapit|db=|gene=tetA|acc=", default_db="megares")
    assert empty.database == "megares"


def test_decode_absent_func_is_empty() -> None:
    """Given a header without a func key, When decoded, Then function is the
    empty string (func is intentionally optional incidental metadata)."""
    header = decode_seqid("gapit|db=ncbi|gene=tetA|acc=NC_1", default_db="ncbi")
    assert header.function == ""


def test_decode_absent_acc_raises() -> None:
    """Given a header without an acc key, When decoded, Then DatabaseError
    HEADER_MALFORMED with machine-stable reason missing_acc."""
    with pytest.raises(DatabaseError) as excinfo:
        decode_seqid("gapit|db=ncbi|gene=tetA", default_db="ncbi")
    assert excinfo.value.code == "HEADER_MALFORMED"
    assert excinfo.value.context["reason"] == "missing_acc"


def test_decode_absent_db_raises() -> None:
    """Given a header without a db key, When decoded, Then DatabaseError
    HEADER_MALFORMED with machine-stable reason missing_db."""
    with pytest.raises(DatabaseError) as excinfo:
        decode_seqid("gapit|gene=tetA|acc=NC_1", default_db="ncbi")
    assert excinfo.value.code == "HEADER_MALFORMED"
    assert excinfo.value.context["reason"] == "missing_db"


def test_decode_empty_acc_value_allowed() -> None:
    """Given a header with an empty acc value (unpublished sequence), When
    decoded, Then accession is the empty string without error."""
    header = decode_seqid("gapit|db=ncbi|gene=tetA|acc=|func=", default_db="ncbi")
    assert header.accession == ""


def test_decode_missing_gene_raises() -> None:
    """Given a header without a gene key, or with an empty gene value, When
    decoded, Then DatabaseError HEADER_MALFORMED."""
    for seqid in ("gapit|db=ncbi", "gapit|db=ncbi|gene="):
        with pytest.raises(DatabaseError) as excinfo:
            decode_seqid(seqid, default_db="ncbi")
        assert excinfo.value.code == "HEADER_MALFORMED"


def test_decode_bad_escape_raises() -> None:
    """Given a value with '%' not followed by two hex digits (%zz, truncated %3),
    When decoded, Then DatabaseError HEADER_MALFORMED."""
    for seqid in ("gapit|gene=a%zz", "gapit|gene=a%3", "gapit|gene=%"):
        with pytest.raises(DatabaseError) as excinfo:
            decode_seqid(seqid, default_db="ncbi")
        assert excinfo.value.code == "HEADER_MALFORMED"


def test_decode_duplicate_key_raises() -> None:
    """Given a header repeating the gene key, When decoded, Then DatabaseError
    HEADER_MALFORMED."""
    with pytest.raises(DatabaseError) as excinfo:
        decode_seqid("gapit|gene=a|gene=b", default_db="ncbi")
    assert excinfo.value.code == "HEADER_MALFORMED"


def test_decode_segment_without_equals_raises() -> None:
    """Given a tagged segment with no '=' (not a key/value), When decoded, Then
    DatabaseError HEADER_MALFORMED."""
    with pytest.raises(DatabaseError) as excinfo:
        decode_seqid("gapit|db=ncbi|gene", default_db="ncbi")
    assert excinfo.value.code == "HEADER_MALFORMED"


def test_decode_error_carries_seqid_context() -> None:
    """Given a malformed header, When decoded, Then the typed error envelope has
    code HEADER_MALFORMED, exit code 4, and context with the seqid."""
    with pytest.raises(DatabaseError) as excinfo:
        decode_seqid("gapit|gene=a%zz", default_db="ncbi")
    assert excinfo.value.code == "HEADER_MALFORMED"
    assert excinfo.value.exit_code == 4
    assert excinfo.value.context["seqid"] == "gapit|gene=a%zz"
    assert excinfo.value.context["reason"]


def test_decode_legacy_tilde_header_delegates() -> None:
    """Given a legacy ~~~ header (full and partial), When decoded, Then the
    result matches db.py parse_db_header exactly."""
    full = "ncbi~~~tetA~~~NC_000913.3~~~TETRACYCLINE"
    partial = "vfdb~~~toxA"
    assert decode_seqid(full, default_db="ncbi") == parse_db_header(full, "ncbi")
    assert decode_seqid(partial, default_db="ncbi") == parse_db_header(partial, "ncbi")


def test_decode_legacy_no_separator_delegates() -> None:
    """Given a seqid with no separator at all, When decoded, Then the whole id
    is the gene under the default database (legacy rule)."""
    header = decode_seqid("tetA(1)", default_db="resfinder")
    assert header.gene == "tetA(1)"
    assert header.database == "resfinder"
    assert header.accession == ""
    assert header.function == ""


def test_is_gapit_header_sniff() -> None:
    """Given various seqids, When sniffed, Then only a literal 'gapit|' prefix
    counts as the new format."""
    assert is_gapit_header("gapit|db=ncbi|gene=tetA")
    assert is_gapit_header("gapit|")
    assert not is_gapit_header("ncbi~~~tetA")
    assert not is_gapit_header("gapit")
    assert not is_gapit_header("xgapit|db=ncbi")
