"""Tests for the plasmidfinder provider transform (Wave B8) — offline only.

The bitbucket HEAD.zip carries an arbitrary top-level directory (the repo
commit hash); the tests assemble that archive from committed ``.fsa`` fixtures
with the stdlib ``zipfile`` — no network, no fetch_provider, no makeblastdb
(the orchestrator gates the assembled pipeline once).

Perl ground truth (abricate-get_db 1.4.0 ``get_plasmidfinder``, verified
against the parity env): ``demo_repB_1_NC_011111.1`` yields gene
``demo_repB_1`` + accession ``NC_011111.1`` — the ``_1`` copy number STAYS on
the gene; only the accession suffix is split off. A non-matching id keeps
itself as the gene with an empty accession. The product is always the
ORIGINAL full id (upstream sets DESC before the regex munging).
"""

import zipfile
from collections.abc import Mapping
from pathlib import Path

from gapit.providers.common import Provider
from gapit.providers.plasmidfinder import PROVIDER, transform
from gapit.records import Record

FIXTURES = Path(__file__).parent / "data" / "providers" / "plasmidfinder"
TOP = "fake-pf-xyz"  # arbitrary bitbucket-style top-level directory

# 60 bp per record (repo fixture convention; nothing here builds an index).
SEQ_REPB = "ACGTTGCAAG" * 6
SEQ_REPA = "TTGCAAGGCC" * 6
SEQ_BARE = "GGCCTTAAGC" * 6

EXPECTED = (
    Record(
        db="plasmidfinder",
        gene="demo_repB_1",
        accession="NC_011111.1",
        function=("replicon",),
        product="demo_repB_1_NC_011111.1",
        sequence=SEQ_REPB,
        source_id="demo_repB_1_NC_011111.1",
    ),
    Record(
        db="plasmidfinder",
        gene="demo_repA_2",
        accession="FAKE123",
        function=("replicon",),
        product="demo_repA_2_FAKE123",
        sequence=SEQ_REPA,
        source_id="demo_repA_2_FAKE123",
    ),
    Record(
        db="plasmidfinder",
        gene="bare_header",
        accession="",
        function=("replicon",),
        product="bare_header",
        sequence=SEQ_BARE,
        source_id="bare_header",
    ),
)


def _build_zip(workdir: Path, members: Mapping[str, str]) -> None:
    """Assemble ``workdir/HEAD.zip`` with the fake top dir and fixed timestamps
    (deterministic bytes; mirrors the bitbucket download layout)."""
    with zipfile.ZipFile(workdir / "HEAD.zip", "w") as archive:
        for name, content in members.items():
            info = zipfile.ZipInfo(f"{TOP}/{name}", date_time=(1980, 1, 1, 0, 0, 0))
            archive.writestr(info, content)


def test_transform_parses_archive_into_exact_records(tmp_path: Path) -> None:
    """Given a HEAD.zip built from the committed fixtures (two .fsa members
    under an arbitrary top-level directory), When transformed, Then the exact
    Records appear in sorted-glob order: matching ids split the accession off
    the gene (the ``_1``/``_2`` copy number stays on the gene), the
    non-matching id keeps itself as gene with an empty accession, and every
    product and source_id is the ORIGINAL full id."""
    _build_zip(
        tmp_path,
        {
            name: (FIXTURES / name).read_text(encoding="utf-8")
            for name in ("plasmids_a.fsa", "plasmids_b.fsa")
        },
    )
    assert tuple(transform(tmp_path)) == EXPECTED


def test_transform_strips_trailing_underscore_runs_from_matched_id(tmp_path: Path) -> None:
    """Given an id whose regex capture ends in an underscore run
    (``demo_repX_9__ABC123`` — the separator consumes one ``_`` of the pair),
    When transformed, Then the gene is ``demo_repX_9``: the captured id with
    trailing underscores stripped (perl ``s/_+$//g``)."""
    _build_zip(tmp_path, {"inline.fsa": f">demo_repX_9__ABC123\n{'CCAAGGTTAC' * 6}\n"})
    (record,) = transform(tmp_path)
    assert record.gene == "demo_repX_9"
    assert record.accession == "ABC123"


def test_provider_metadata_matches_the_frozen_contract() -> None:
    """Given the plasmidfinder provider module, When inspected, Then it exposes
    the B0 Provider contract with the pinned name, description, source URL,
    dbtype, and transform wiring."""
    assert isinstance(PROVIDER, Provider)
    assert PROVIDER.name == "plasmidfinder"
    assert PROVIDER.description == "CGE PlasmidFinder replicons"
    assert PROVIDER.source_urls == (
        "https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/HEAD.zip",
    )
    assert PROVIDER.dbtype == "nucl"
    assert PROVIDER.transform is transform
