"""Tests for the ecoh provider transform (Wave B10) — offline, fixture-only.

Upstream headers follow the srst2 convention
``>cluster__gene__allele__seq_id acc;piece;piece`` (abricate-get_db
``get_ecoh``): the id splits on ``__`` into 4 parts with the gene taken from
part 3 (the allele, e.g. ``fliC-H1``); the description splits on ``;`` into
an accession plus a space-joined product. Ids lacking 4 ``__``-parts are
skipped. These tests never touch the network, fetch_provider, or
makeblastdb — the orchestrator gates the assembled pipeline once.
"""

from pathlib import Path

from gapit.providers.common import Provider
from gapit.providers.ecoh import PROVIDER, transform
from gapit.records import Record

WORKDIR = Path(__file__).parent / "data" / "providers" / "ecoh"

# 60 bp per record (repo fixture convention; nothing here builds an index).
SEQ_FLIC = "ACGTTGCAAG" * 6
SEQ_WZX = "TTGCAAGGCC" * 6
SEQ_SYN = "AACCGGTTAA" * 6

EXPECTED = (
    Record(
        db="ecoh",
        gene="fliC-H1",
        accession="FAKE0001.1",
        function=("H-antigen",),
        product="flagellin H1",
        sequence=SEQ_FLIC,
        source_id="1__fliC__fliC-H1__1",
    ),
    Record(
        db="ecoh",
        gene="wzx-O41",
        accession="FAKE0002.1",
        function=("O-antigen",),
        product="O antigen flippase O41",
        sequence=SEQ_WZX,
        source_id="8__wzx__wzx-O41__246",
    ),
    Record(
        db="ecoh",
        gene="syn-ANTIGEN",
        accession="FAKE0004.1",
        function=("antigen",),
        product="clearly synthetic fallback antigen marker",
        sequence=SEQ_SYN,
        source_id="30__syn__syn-ANTIGEN__7",
    ),
)


def test_transform_parses_fixture_into_exact_records() -> None:
    """Given the committed EcOH fixture (one record with a two-piece
    description, one with three pieces, one clearly-synthetic fallback), When
    transformed, Then the exact Records are yielded: db ecoh, gene from id
    part 3, accession from description piece 1, product as the space-joined
    remainder, sequence and source_id carried verbatim, function derived
    from the allele prefix."""
    assert tuple(transform(WORKDIR)) == EXPECTED


def test_function_category_follows_the_allele_prefix() -> None:
    """Given fixture alleles covering every derivation branch, When
    transformed, Then function is allele-prefix derived: fliC* → H-antigen,
    wzx/wzy/wzt/wzm* → O-antigen, anything else (the synthetic syn-ANTIGEN)
    → the "antigen" fallback."""
    records = {record.gene: record for record in transform(WORKDIR)}
    assert records["fliC-H1"].function == ("H-antigen",)
    assert records["wzx-O41"].function == ("O-antigen",)
    assert records["syn-ANTIGEN"].function == ("antigen",)


def test_transform_skips_ids_without_four_double_underscore_parts() -> None:
    """Given a fixture record whose id contains no '__' separator at all,
    When transformed, Then no Record is produced for it — only the three
    well-formed records survive."""
    records = tuple(transform(WORKDIR))
    assert len(records) == 3
    assert all(record.source_id != "malformed-id" for record in records)


def test_provider_metadata_matches_the_frozen_contract() -> None:
    """Given the ecoh provider module, When inspected, Then it exposes the
    B0 Provider contract with the pinned name, description, source URL,
    dbtype, and transform wiring."""
    assert isinstance(PROVIDER, Provider)
    assert PROVIDER.name == "ecoh"
    assert PROVIDER.description == "E. coli O and H antigens (srst2 EcOH)"
    assert PROVIDER.source_urls == (
        "https://raw.githubusercontent.com/katholt/srst2/master/data/EcOH.fasta",
    )
    assert PROVIDER.dbtype == "nucl"
    assert PROVIDER.transform is transform
