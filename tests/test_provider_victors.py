"""Tests for the victors provider (Wave B4) — transform only, offline.

Fixtures under tests/data/providers/victors mirror the two phidias.us
downloads exactly as they land in a fetch workdir (.php URL basenames
carrying .ffn/.faa content). No network, no fetch_provider, no build —
the B0 pipeline is already covered by test_providers_common.py.
"""

from pathlib import Path

from gapit.providers.victors import PROVIDER, transform
from gapit.records import Record

DATA = Path(__file__).parent / "data" / "providers" / "victors"

SEQ_MAP_BRACKET = "ACGTTGCA" * 7 + "ACGT"  # 60 bp, minimap2-indexable
SEQ_MAP_CLEAN = "GATCCGTA" * 7 + "GATC"
SEQ_FALLBACK = "TTACGGCA" * 7 + "TTAC"

# gi 115534244: map hit, .faa product carries a strain bracket.
EXPECTED_BRACKET = Record(
    db="victors",
    gene="gi|115534244:5211-5733",
    sequence=SEQ_MAP_BRACKET,
    accession="YP_783826.1",
    function=("virulence",),
    product="hypothetical protein pCJ01p4",
    source_id="115534244",
)
# gi 115534245: map hit, clean .faa product without a bracket.
EXPECTED_CLEAN = Record(
    db="victors",
    gene="gi|115534245:5734-6255",
    sequence=SEQ_MAP_CLEAN,
    accession="YP_783827.1",
    function=("virulence",),
    product="conjugal transfer protein TraG",
    source_id="115534245",
)
# gi 115534241: absent from the .faa -> coordinates fallback.
EXPECTED_FALLBACK = Record(
    db="victors",
    gene="gi|115534241:2616-3152",
    sequence=SEQ_FALLBACK,
    accession="gi|115534241:2616-3152",
    function=("virulence",),
    product="hypothetical protein",
    source_id="115534241",
)


def test_provider_metadata_matches_spec() -> None:
    """Given the victors PROVIDER, When inspected, Then name, description,
    both source URLs in download order, and dbtype match the Wave B4 spec
    (frozen B0 contract)."""
    assert PROVIDER.name == "victors"
    assert PROVIDER.description == "Victors virulence factors"
    assert PROVIDER.source_urls == (
        "http://phidias.us/victors/downloads/gen_downloads.php",
        "http://phidias.us/victors/downloads/gen_downloads_protein.php",
    )
    assert PROVIDER.dbtype == "nucl"


def test_transform_yields_exact_records() -> None:
    """Given the fixture downloads (3 protein headers, 3 CDS records, one
    protein-only gi), When transformed, Then exactly the three expected
    Records come out, in .ffn order, all with db victors."""
    assert tuple(transform(DATA)) == (
        EXPECTED_BRACKET,
        EXPECTED_CLEAN,
        EXPECTED_FALLBACK,
    )


def test_transform_maps_faa_accession_and_product() -> None:
    """Given an .ffn record whose gi has a clean .faa header, When
    transformed, Then both accession and product come from the .faa map."""
    by_source = {record.source_id: record for record in transform(DATA)}
    assert by_source["115534245"].accession == "YP_783827.1"
    assert by_source["115534245"].product == "conjugal transfer protein TraG"


def test_transform_truncates_product_at_strain_bracket() -> None:
    """Given a .faa header whose product ends in a `` [strain]`` bracket,
    When transformed, Then the product stops before the bracket (upstream
    ``([^[]+)`` semantics)."""
    by_source = {record.source_id: record for record in transform(DATA)}
    assert by_source["115534244"].product == "hypothetical protein pCJ01p4"


def test_transform_falls_back_when_gi_missing_from_faa() -> None:
    """Given an .ffn record whose gi has NO .faa header, When transformed,
    Then accession falls back to the gi|N:start-stop form and product to
    'hypothetical protein'."""
    fallback = {record.source_id: record for record in transform(DATA)}["115534241"]
    assert fallback.accession == "gi|115534241:2616-3152"
    assert fallback.product == "hypothetical protein"


def test_transform_skips_records_without_gi_coordinates(tmp_path: Path) -> None:
    """Given an .ffn holding a record whose id carries no gi|N:s-e part,
    When transformed, Then that record is skipped (upstream would build it
    from the PREVIOUS match's stale capture variables) and the well-formed
    record survives."""
    (tmp_path / "gen_downloads_protein.php").write_text(
        ">gi|115534245|ref|YP_783827.1| conjugal transfer protein TraG\nMKV\n",
        encoding="utf-8",
    )
    (tmp_path / "gen_downloads.php").write_text(
        f">gi|115534245:5734-6255 CDS pCJ01p05\n{SEQ_MAP_CLEAN}\n"
        f">plasmid_knob no gi coordinates here\n{SEQ_FALLBACK}\n",
        encoding="utf-8",
    )
    records = tuple(transform(tmp_path))
    assert [record.gene for record in records] == ["gi|115534245:5734-6255"]
