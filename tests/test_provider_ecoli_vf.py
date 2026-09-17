"""Tests for the ecoli_vf provider (Wave B6) — transform only, offline.

The upstream ``repaired_ecoli_vfs_shortnames.ffn`` is simulated by a
committed plain fixture (the download basename, uncompressed): no network,
no fetch_provider, no build — the B0 pipeline is already covered by
test_providers_common.py.
"""

from pathlib import Path

import pytest

from gapit.providers.ecoli_vf import PROVIDER, transform
from gapit.records import Record

FIXTURE = Path(__file__).parent / "data" / "providers" / "ecoli_vf"
SOURCE_FILE = "repaired_ecoli_vfs_shortnames.ffn"
SOURCE_URL = (
    "https://github.com/phac-nml/ecoli_vf/raw/master/data/repaired_ecoli_vfs_shortnames.ffn"
)
SEQ_ESPX = "ACGT" * 15  # 60 bp, minimap2-indexable per Wave A3
SEQ_CVII = "CAGT" * 15
SEQ_PLAIN = "GTCA" * 15


@pytest.fixture
def workdir(tmp_path: Path) -> Path:
    """A download workdir holding the fixture under the name the provider's
    source URL would land it as (plain file — byte copy, no decompression)."""
    (tmp_path / SOURCE_FILE).write_bytes((FIXTURE / SOURCE_FILE).read_bytes())
    return tmp_path


@pytest.fixture
def records(workdir: Path) -> tuple[Record, ...]:
    """All Records the fixture yields, in file order."""
    return tuple(transform(workdir))


def test_provider_metadata_matches_spec() -> None:
    """Given the ecoli_vf PROVIDER, When inspected, Then name, description,
    source URLs, and dbtype match the Wave B6 spec (frozen B0 contract)."""
    assert PROVIDER.name == "ecoli_vf"
    assert PROVIDER.description == "E. coli virulence factors (phac-nml)"
    assert PROVIDER.source_urls == (SOURCE_URL,)
    assert PROVIDER.dbtype == "nucl"


def test_id_paren_desc_override_and_two_trailing_brackets(records: tuple[Record, ...]) -> None:
    """Given a record with an accession paren in the id, a desc-paren gene
    override, and TWO trailing bracket groups, When transformed, Then both
    brackets are stripped in one pass, the override renames the gene, the
    paren content becomes the accession, and source_id keeps the base id."""
    assert records[0] == Record(
        db="ecoli_vf",
        gene="espX",
        sequence=SEQ_ESPX,
        accession="gi:999123",
        function=("virulence",),
        product="Demo protein",
        source_id="DEMO000748",
    )


def test_missing_id_paren_falls_back_to_base_accession(records: tuple[Record, ...]) -> None:
    """Given a record whose id has NO paren suffix, When transformed, Then
    the accession falls back to the base id ($2 || $1) while the desc-paren
    still renames the gene and the trailing bracket group is stripped."""
    assert records[1] == Record(
        db="ecoli_vf",
        gene="cvii",
        sequence=SEQ_CVII,
        accession="SPG000142",
        function=("virulence",),
        product="Fake coli cvi cvaC operon.",
        source_id="SPG000142",
    )


def test_plain_record_keeps_base_gene(records: tuple[Record, ...]) -> None:
    """Given a record with no parens or brackets anywhere, When transformed,
    Then gene and accession are both the base id and the description is the
    product verbatim."""
    assert records[2] == Record(
        db="ecoli_vf",
        gene="PLAIN0001",
        sequence=SEQ_PLAIN,
        accession="PLAIN0001",
        function=("virulence",),
        product="hypothetical demo protein",
        source_id="PLAIN0001",
    )


def test_transform_skips_records_whose_id_fails_the_pattern(tmp_path: Path) -> None:
    """Given the fixture plus a record whose id cannot match
    ``^(\\w+)(?:\\((.*?)\\))?$``, When transformed, Then only the three
    well-formed records survive (upstream die()s on such ids; we skip —
    documented deviation)."""
    (tmp_path / SOURCE_FILE).write_bytes((FIXTURE / SOURCE_FILE).read_bytes())
    with (tmp_path / SOURCE_FILE).open("a", encoding="utf-8") as handle:
        handle.write(f">BAD-RECORD unparseable id keeps a hyphen\n{'AAAA' * 15}\n")
    survivors = tuple(transform(tmp_path))
    assert [record.gene for record in survivors] == ["espX", "cvii", "PLAIN0001"]


def test_empty_id_paren_falls_back_to_base_accession(tmp_path: Path) -> None:
    """Given a record whose id paren group is EMPTY, When transformed, Then
    the accession falls back to the base id — perl's ``$2 || $1`` is
    falsy-based, so an empty capture does not become the accession."""
    (tmp_path / SOURCE_FILE).write_text(
        f">EMPTY() product line only\n{SEQ_ESPX}\n", encoding="utf-8"
    )
    assert tuple(transform(tmp_path)) == (
        Record(
            db="ecoli_vf",
            gene="EMPTY",
            sequence=SEQ_ESPX,
            accession="EMPTY",
            function=("virulence",),
            product="product line only",
            source_id="EMPTY",
        ),
    )
