"""Tests for the resfinder provider transform (Wave B3) — offline only.

The bitbucket HEAD.zip carries an arbitrary top-level directory; the tests
wrap the committed fixtures into that archive with the stdlib zipfile — no
network, no fetch_provider, no makeblastdb (the orchestrator gates the
assembled pipeline once for all providers).

Perl ground truth (abricate-get_db 1.4.0 ``get_resfinder``):
``phenotypes.txt`` col2 supplies the function categories (split on
comma+whitespace, unknown/notes/none markers filtered, LAST row per gene
wins); ``.fsa`` ids split into gene prefix, copy number, and accession
(``demoA_1_FAKE0001`` -> gene ``demoA_1``, accession ``FAKE0001``, product
``demoA``); glued headers (abricate issue #62) are repaired before parsing.
An id with no ``_<digits>_<accession>`` structure is skipped — upstream
would emit it with undef fields.
"""

import zipfile
from pathlib import Path

import pytest

from gapit.errors import DatabaseError
from gapit.providers.common import Provider
from gapit.providers.resfinder import PROVIDER, transform
from gapit.records import Record

FIXTURES = Path(__file__).parent / "data" / "providers" / "resfinder"
TOP = "fake-resfinder-abc123"  # arbitrary bitbucket-style top-level directory

ISSUE62_TAIL = "GCTTTAAATTGGAAAAAAGATAGTCAAAACTCTTTAA"  # glued before >demoE
# 60 bp per record (repo fixture convention; nothing here builds an index).
SEQ_DEMOB = "TTTAAACCCGGG" * 5
SEQ_DEMOC = "GGGCAT" * 10
SEQ_DEMOA = "ACGTT" * 12 + ISSUE62_TAIL  # the glued line ends demoA's sequence
SEQ_DEMOE = "CCCATTCGGA" * 6

EXPECTED = (
    Record(  # aminoglycoside.fsa: class "unknown" filters to empty function
        db="resfinder",
        gene="demoB_1",
        accession="FAKE0002",
        function=(),
        product="demoB",
        sequence=SEQ_DEMOB,
        source_id="demoB_1_FAKE0002",
    ),
    Record(  # aminoglycoside.fsa: LAST phenotypes row for demoC wins
        db="resfinder",
        gene="demoC_2",
        accession="FAKE0004",
        function=("Class LAST",),
        product="demoC",
        sequence=SEQ_DEMOC,
        source_id="demoC_2_FAKE0004",
    ),
    Record(  # beta-lactam.fsa: multi-class row splits on comma-space
        db="resfinder",
        gene="demoA_1",
        accession="FAKE0001",
        function=("Class A", "Class B"),
        product="demoA",
        sequence=SEQ_DEMOA,
        source_id="demoA_1_FAKE0001",
    ),
    Record(  # beta-lactam.fsa: the repaired glued header is its own record
        db="resfinder",
        gene="demoE_1",
        accession="FAKE0005",
        function=(),  # demoE has no phenotypes.txt row
        product="demoE",
        sequence=SEQ_DEMOE,
        source_id="demoE_1_FAKE0005",
    ),
)


def _head_zip(workdir: Path) -> Path:
    """Assemble ``workdir/HEAD.zip`` from the committed fixtures under the
    fake top-level directory (the bitbucket download layout)."""
    archive_path = workdir / "HEAD.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        for name in ("phenotypes.txt", "aminoglycoside.fsa", "beta-lactam.fsa"):
            archive.write(FIXTURES / name, arcname=f"{TOP}/{name}")
    return archive_path


def test_transform_parses_archive_into_exact_records(tmp_path: Path) -> None:
    """Given a HEAD.zip wrapping phenotypes.txt and two .fsa fixtures under an
    arbitrary top-level directory, When transformed, Then the exact Records
    appear in sorted-file order: gene is prefix_copy, accession is the id
    suffix, product is the bare prefix, and function comes from the gene's
    phenotypes row (or is empty without one)."""
    _head_zip(tmp_path)
    assert tuple(transform(tmp_path)) == EXPECTED


def test_transform_repairs_glued_inline_header(tmp_path: Path) -> None:
    """Given beta-lactam.fsa glues demoE's header onto the tail of demoA's
    final sequence line (abricate issue #62), When transformed, Then demoE
    parses as its OWN record with a clean sequence, and no '>' survives
    inside any sequence."""
    _head_zip(tmp_path)
    records = {record.gene: record for record in transform(tmp_path)}
    assert records["demoE_1"].sequence == SEQ_DEMOE
    assert records["demoE_1"].source_id == "demoE_1_FAKE0005"
    assert records["demoA_1"].sequence == SEQ_DEMOA  # tail, no '>' char
    assert all(">" not in record.sequence for record in records.values())


def test_transform_keeps_last_phenotypes_row_per_gene(tmp_path: Path) -> None:
    """Given phenotypes.txt lists demoC twice (Class OLD, then Class LAST),
    When transformed, Then the demoC record carries only the LAST row's
    classes (upstream's plain hash overwrite)."""
    _head_zip(tmp_path)
    (democ,) = [record for record in transform(tmp_path) if record.gene == "demoC_2"]
    assert democ.function == ("Class LAST",)


def test_transform_filters_unknown_class_to_empty_function(tmp_path: Path) -> None:
    """Given demoB's phenotypes row has the single class 'unknown', When
    transformed, Then every piece is filtered out and the record's
    function tuple is empty (not the literal string 'unknown')."""
    _head_zip(tmp_path)
    (demob,) = [record for record in transform(tmp_path) if record.gene == "demoB_1"]
    assert demob.function == ()


def test_transform_skips_ids_without_copy_and_accession(tmp_path: Path) -> None:
    """Given a fixture record id with no ``_<digits>_<accession>`` structure
    ('no_copy_number'), When transformed, Then no Record is emitted for it —
    the documented deviation from upstream's undef-field record."""
    _head_zip(tmp_path)
    records = tuple(transform(tmp_path))
    assert len(records) == 4
    assert all(record.source_id != "no_copy_number" for record in records)


def test_transform_raises_typed_error_without_phenotypes(tmp_path: Path) -> None:
    """Given a HEAD.zip containing .fsa members but no phenotypes.txt, When
    transformed, Then a typed PROVIDER_INVALID DatabaseError is raised
    (upstream dies reading the missing file rather than guessing classes)."""
    with zipfile.ZipFile(tmp_path / "HEAD.zip", "w") as archive:
        archive.write(FIXTURES / "aminoglycoside.fsa", arcname=f"{TOP}/aminoglycoside.fsa")
    with pytest.raises(DatabaseError) as raised:
        tuple(transform(tmp_path))
    assert raised.value.code == "PROVIDER_INVALID"


def test_provider_metadata_matches_the_frozen_contract() -> None:
    """Given the resfinder provider module, When inspected, Then it exposes
    the B0 Provider contract with the pinned name, description, source URL,
    dbtype, and transform wiring."""
    assert isinstance(PROVIDER, Provider)
    assert PROVIDER.name == "resfinder"
    assert PROVIDER.description == "CGE ResFinder acquired resistance genes"
    assert PROVIDER.source_urls == (
        "https://bitbucket.org/genomicepidemiology/resfinder_db/get/HEAD.zip",
    )
    assert PROVIDER.dbtype == "nucl"
    assert PROVIDER.transform is transform
