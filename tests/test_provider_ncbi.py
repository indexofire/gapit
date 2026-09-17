"""Tests for the NCBI AMRFinderPlus provider (providers/ncbi.py, Wave B1).

Pure transform tests over synthetic fixtures under tests/data/providers/ncbi/:
no network, no fetch_provider, no build (Wave B0/C concerns). Every
abricate-get_db ``get_ncbi`` rule is pinned: the 7-field ``|`` id split, the
fp/fn fusion filter, the ``.1`` accession append for unversioned accessions,
the core/AMR/AMR catalog gating, subclass split into function categories, and
the underscore->space product rewrite.
"""

from pathlib import Path

import pytest

from gapit.providers.ncbi import PROVIDER, transform
from gapit.records import Record

WORKDIR = Path(__file__).parent / "data" / "providers" / "ncbi"

SEQ_PLAIN = "ACGT" * 15
SEQ_FUSION = "ACCC" * 15
SEQ_APPEND = "AGGG" * 15
SEQ_SCOPE = "CATG" * 15
SEQ_TYPE = "TGCA" * 15
SEQ_MISSING = "GGCC" * 15
SEQ_MULTI = "TTAA" * 15

# Expected output of transform(WORKDIR) in fixture FASTA order. The append
# keeper's catalog key is the versioned accession only (fasta acc NG_900003),
# and its accession comes from the TSV row, not the fasta header.
EXPECTED_RECORDS = (
    Record(
        db="ncbi",
        gene="fakeGenePlain",
        sequence=SEQ_PLAIN,
        accession="NG_900001.1",
        function=("fakeclass",),
        product="fake product alpha",
        source_id="WP_900000001.1",
    ),
    Record(
        db="ncbi",
        gene="fakeGeneAppend",
        sequence=SEQ_APPEND,
        accession="NG_900003.1",
        function=("solo",),
        product="fake product gamma",
        source_id="WP_900000003.1",
    ),
    Record(
        db="ncbi",
        gene="fakeGeneMulti",
        sequence=SEQ_MULTI,
        accession="NG_900007.1",
        function=("a", "b"),
        product="fake product eta",
        source_id="WP_900000007.1",
    ),
)

SKIPPED_GENES = (
    "fakeGeneFusion",
    "fakeGeneScope",
    "fakeGeneType",
    "fakeGeneMissing",
)


def test_transform_yields_exactly_the_qualifying_records() -> None:
    """Given the synthetic AMR_CDS.fa + ReferenceGeneCatalog.txt fixtures,
    When transform runs over the download workdir, Then exactly the three
    qualifying Records are yielded in FASTA order with every field mapped
    the get_ncbi way: db fixed to ncbi, accession taken from the catalog row
    (proving the '.1' append for NG_900003), subclass split on '/' into
    function categories, product underscores rewritten to spaces, and
    source_id from the id's first field."""
    assert tuple(transform(WORKDIR)) == EXPECTED_RECORDS


@pytest.mark.parametrize("gene", SKIPPED_GENES)
def test_transform_skips_record_when_one_rule_excludes_it(gene: str) -> None:
    """Given a fixture FASTA record excluded by exactly one get_ncbi rule
    (fp=0 fusion, catalog scope 'bonus', catalog type STRESS, accession with
    no catalog row), When transform runs, Then its unique gene is absent
    from every yielded record — passively skipped, no error, like upstream."""
    assert gene not in {record.gene for record in transform(WORKDIR)}


def test_provider_metadata_matches_abricate_get_db() -> None:
    """Given the frozen Wave B0 Provider contract, When the ncbi module is
    imported, Then PROVIDER carries the exact get_ncbi name, description,
    AMRFinderPlus latest URL pair, and nucl dbtype, wired to this module's
    transform."""
    base = (
        "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/"
        "AMRFinderPlus/database/latest"
    )
    assert PROVIDER.name == "ncbi"
    assert PROVIDER.description == "NCBI AMRFinderPlus (reference finder) curated AMR"
    assert PROVIDER.source_urls == (f"{base}/AMR_CDS.fa", f"{base}/ReferenceGeneCatalog.txt")
    assert PROVIDER.dbtype == "nucl"
    assert PROVIDER.transform is transform
