"""Tests for the megares provider transform (Wave B9) — offline, fixture-only.

MEGARes v3 headers are pipe-separated with up to 6 fields
``id|type|class|mech|group|note`` (abricate-get_db ``get_megares``): the note
field exists only on records requiring SNP confirmation, which are skipped
(SPEC §8); keepers carry 5 non-empty fields. The download archive is built
in-test from the committed plain fasta (no binary fixtures in git), nested
one directory deep like the real zip. These tests never touch the network,
fetch_provider, or makeblastdb — the orchestrator gates the assembled
pipeline once.
"""

from pathlib import Path
from zipfile import ZipFile

from gapit.providers.common import Provider
from gapit.providers.megares import PROVIDER, transform
from gapit.records import Record

FIXTURE = Path(__file__).parent / "data" / "providers" / "megares" / "megares_drugs_demo.fasta"
ARCHIVE = "megares_v3.00.zip"

# 60 bp per record (repo fixture convention; nothing here builds an index).
SEQ_DEMO1 = "ACGTTGCAAG" * 6

EXPECTED = (
    Record(
        db="megares",
        gene="DEMO1",
        accession="MEG_9100",
        function=("demo_class",),
        product="Drugs:demo_class:demo_mech:DEMO1",
        sequence=SEQ_DEMO1,
        source_id="MEG_9100",
    ),
)


def _workdir(tmp_path: Path) -> Path:
    """A download workdir holding megares_v3.00.zip with the fixture fasta
    nested one directory deep (the real archive layout; upstream's
    ``unzip -j`` flattens, the transform globs recursively)."""
    workdir = tmp_path / "work"
    workdir.mkdir()
    with ZipFile(workdir / ARCHIVE, "w") as archive:
        archive.write(FIXTURE, "megares_v3.00/megares_drugs_demo.fasta")
    return workdir


def test_transform_extracts_archive_and_parses_exact_records(tmp_path: Path) -> None:
    """Given a workdir with megares_v3.00.zip containing the demo fasta,
    When transformed, Then exactly the keeper Record is yielded: db megares,
    gene from field 5, accession and source_id from field 1, product as the
    colon-join of fields 2-5, function carrying the class field (field 3) as
    a 1-tuple, sequence verbatim."""
    assert tuple(transform(_workdir(tmp_path))) == EXPECTED


def test_transform_skips_snp_confirmation_and_short_ids(tmp_path: Path) -> None:
    """Given the demo fasta's SNP-confirmation record (non-empty 6th field)
    and its 5-field record with an empty group, When transformed, Then
    neither appears in the output: one record total, no DEMO2 gene, no
    MEG_9101 or MEG_9102 source ids."""
    records = tuple(transform(_workdir(tmp_path)))
    assert len(records) == 1
    assert all(record.gene != "DEMO2" for record in records)
    assert all(record.source_id not in {"MEG_9101", "MEG_9102"} for record in records)


def test_provider_metadata_matches_the_frozen_contract() -> None:
    """Given the megares provider module, When inspected, Then it exposes the
    B0 Provider contract with the pinned name, description, source URL,
    dbtype, and transform wiring."""
    assert isinstance(PROVIDER, Provider)
    assert PROVIDER.name == "megares"
    assert PROVIDER.description == "MEGARes antimicrobial resistance genes"
    assert PROVIDER.source_urls == ("https://www.meglab.org/downloads/megares_v3.00.zip",)
    assert PROVIDER.dbtype == "nucl"
    assert PROVIDER.transform is transform
