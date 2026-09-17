"""Tests for the vfdb provider (Wave B5) — transform only, offline.

The upstream ``VFDB_setA_nt.fas.gz`` is simulated by gzip-compressing a
committed plain fixture at test time (gzip stdlib): no network, no
fetch_provider, no build — the B0 pipeline is already covered by
test_providers_common.py.
"""

import gzip
from pathlib import Path

import pytest

from gapit.providers.vfdb import PROVIDER, transform
from gapit.records import Record

FIXTURE = Path(__file__).parent / "data" / "providers" / "vfdb" / "VFDB_setA_nt.fas"
SEQ_DEMO = "ACGT" * 15  # 60 bp, minimap2-indexable per Wave A3

EXPECTED = Record(
    db="vfdb",
    gene="demoT",
    sequence=SEQ_DEMO,
    accession="FAKE0001",  # version suffix .2 consumed by the optional group
    function=("virulence",),
    product="(demoT) demo toxin precursor [Demo VF (VF9999)] [Fakeus demos]",
    source_id="DEMO000123(gb|FAKE0001.2)",
)


@pytest.fixture
def workdir(tmp_path: Path) -> Path:
    """A download workdir holding the fixture gzipped under the name the
    provider's source URL would land it as."""
    (tmp_path / "VFDB_setA_nt.fas.gz").write_bytes(gzip.compress(FIXTURE.read_bytes(), mtime=0))
    return tmp_path


def test_provider_metadata_matches_spec() -> None:
    """Given the vfdb PROVIDER, When inspected, Then name, description,
    source URLs, and dbtype match the Wave B5 spec (frozen B0 contract)."""
    assert PROVIDER.name == "vfdb"
    assert PROVIDER.description == "VFDB virulence factors (set A, nucleotide)"
    assert PROVIDER.source_urls == ("http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz",)
    assert PROVIDER.dbtype == "nucl"


def test_transform_yields_exact_record(workdir: Path) -> None:
    """Given the gzipped fixture with one well-formed record, When
    transformed, Then exactly that Record survives: gene renamed to the
    desc-paren name, accession taken from the id parens with the version
    stripped, product the full description, source_id the original id."""
    assert tuple(transform(workdir)) == (EXPECTED,)


def test_transform_skips_records_failing_either_regex(workdir: Path) -> None:
    """Given the fixture also holds a record whose description lacks the
    leading paren group and one whose id lacks the (db|acc) group, When
    transformed, Then both are skipped (upstream would carry stale $1/$2
    into them; we drop the record instead)."""
    records = tuple(transform(workdir))
    assert [record.gene for record in records] == ["demoT"]
    assert all(record.source_id != "DEMO000456(gb|FAKE0002)" for record in records)
    assert all(record.source_id != "DEMO000789" for record in records)


def test_transform_decodes_latin1_bytes_without_error(tmp_path: Path) -> None:
    """Given a gz whose decompressed BYTES carry 0xA0 (latin-1 nbsp) inside
    one record's description AND inside another record's sequence — the
    real-network Wave E failure shape — When transformed, Then no
    UnicodeDecodeError escapes and both Records carry the \xa0 verbatim
    (sequence junk is fetch_provider's N-normalization job, not the
    transform's; pinned here as raw)."""
    seq_desc = b">DEMO000123(gb|FAKE0001.2) (demoT) demo toxin\xa0precursor [Demo VF (VF9999)]\n"
    seq_seq = b">DEMO000456(gb|FAKE0002) (demoU) unicode-free sibling\n"
    latin1_gz = gzip.compress(
        seq_desc + b"ACGT" * 15 + b"\n" + seq_seq + b"ACGT" * 7 + b"\xa0ACGT" * 7 + b"\n",
        mtime=0,
    )
    tmp_path.joinpath("VFDB_setA_nt.fas.gz").write_bytes(latin1_gz)
    records = tuple(transform(tmp_path))
    assert records == (
        Record(
            db="vfdb",
            gene="demoT",
            sequence="ACGT" * 15,
            accession="FAKE0001",
            function=("virulence",),
            product="(demoT) demo toxin\xa0precursor [Demo VF (VF9999)]",
            source_id="DEMO000123(gb|FAKE0001.2)",
        ),
        Record(
            db="vfdb",
            gene="demoU",
            sequence="ACGT" * 7 + "\xa0ACGT" * 7,
            accession="FAKE0002",
            function=("virulence",),
            product="(demoU) unicode-free sibling",
            source_id="DEMO000456(gb|FAKE0002)",
        ),
    )
