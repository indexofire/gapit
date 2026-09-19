"""Tests for the bacmet2 provider (Wave B11) — transform only, offline.

Upstream ``BacMet2_EXP_database.fasta`` is a PROTEIN file — the one
prot-dbtype provider (screening uses blastx) — simulated
by a committed plain fixture: no network, no fetch_provider, no build (the
B0 pipeline is already covered by test_providers_common.py).
"""

from pathlib import Path

import pytest

from gapit.providers.bacmet2 import PROVIDER, transform
from gapit.records import Record

FIXTURE = Path(__file__).parent / "data" / "providers" / "bacmet2" / "BacMet2_EXP_database.fasta"

# Clearly synthetic protein letters (60 aa each, distinct per record).
SEQ_DEMO_C = "MKTAYIAKQRQISFVKDNGTMSFLKSHFSREVLKSFDENLMNEMKRRVAKQRSTQMFVAY"
SEQ_DEMO_F = "MENFNKFSKHQSDLLFEQTKSLVWNTVTNDDLEKVLKDFFTASSSSSTSTSTTQQPPAAW"

EXPECTED_NORMAL = Record(
    db="bacmet2",
    gene="demoC-BAC0098",  # upstream x[1]-x[0] reversal — the point
    sequence=SEQ_DEMO_C,
    accession="sp:P0A502",  # upstream x[2]:x[3] join
    function=("biocide",),
    product="Demo manganese exporting protein",
    source_id="BAC0098|demoC|sp|P0A502",
)

EXPECTED_FALLBACK = Record(
    db="bacmet2",
    gene="demoF-BAC0123",
    sequence=SEQ_DEMO_F,
    accession="tr:Q9X123",
    function=("biocide",),
    product="demoF-BAC0123",  # empty description -> gene (save_fasta DESC || ID)
    source_id="BAC0123|demoF|tr|Q9X123",
)


@pytest.fixture
def workdir(tmp_path: Path) -> Path:
    """A download workdir holding the fixture under the name the provider's
    source URL would land it as."""
    (tmp_path / "BacMet2_EXP_database.fasta").write_bytes(FIXTURE.read_bytes())
    return tmp_path


def test_provider_metadata_matches_spec() -> None:
    """Given the bacmet2 PROVIDER, When inspected, Then name, description,
    source URLs, and dbtype match the Wave B11 spec — dbtype "prot", the
    one protein database in the set."""
    assert PROVIDER.name == "bacmet2"
    assert (
        PROVIDER.description
        == "BacMet2 experimentally confirmed biocide/resistance genes (protein)"
    )
    assert PROVIDER.source_urls == (
        "http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta",
    )
    assert PROVIDER.dbtype == "prot"


def test_transform_yields_exact_records(workdir: Path) -> None:
    """Given the fixture's two well-formed records, When transformed, Then
    exactly those Records survive in file order: gene is the x[1]-x[0]
    reversal, accession the x[2]:x[3] join, product the description or the
    gene when the description is empty, source_id the original id token,
    sequence the raw protein letters."""
    assert tuple(transform(workdir)) == (EXPECTED_NORMAL, EXPECTED_FALLBACK)


def test_transform_skips_records_with_short_ids(workdir: Path) -> None:
    """Given the fixture also holds a record whose id has only 3
    pipe-separated fields, When transformed, Then it is skipped entirely
    (upstream perl would emit undef-indexed garbage; we drop the record)."""
    records = tuple(transform(workdir))
    assert len(records) == 2
    assert all(record.source_id != "BAC0456|shorty|sp" for record in records)


def test_transform_keeps_records_with_extra_id_fields(tmp_path: Path) -> None:
    """Given a record whose id carries the real-BacMet2 5th field (the
    UniProt entry name, e.g. CTPC_MYCTU — see the upstream perl comment),
    When transformed, Then the record is kept: gene and accession come from
    x[0..3] and the extra field is ignored, exactly like upstream."""
    (tmp_path / "BacMet2_EXP_database.fasta").write_text(
        ">BAC0098|ctpC|sp|P0A502|CTPC_MYCTU Probable manganese/zinc-exporting protein\n"
        "MKTAYIAKQRQISFVKDNGTMSFLKSHFSREVLKSFDENLMNEMKRRVAKQRSTQMFVAY\n",
        encoding="utf-8",
    )
    assert tuple(transform(tmp_path)) == (
        Record(
            db="bacmet2",
            gene="ctpC-BAC0098",
            sequence=SEQ_DEMO_C,
            accession="sp:P0A502",
            function=("biocide",),
            product="Probable manganese/zinc-exporting protein",
            source_id="BAC0098|ctpC|sp|P0A502|CTPC_MYCTU",
        ),
    )
