"""Tests for the argannot provider (Wave B7) — transform only, offline.

The real ``ARG-ANNOT_NT_V6_July2019.txt`` ships with stray backslashes that
upstream repairs by stripping them from the raw text before parsing; the
committed fixture reproduces that corruption so the strip is proven. No
network, no fetch_provider, no build (the B0 pipeline is already covered by
test_providers_common.py).
"""

from pathlib import Path

import pytest

from gapit.providers.argannot import PROVIDER, transform
from gapit.records import Record

FIXTURE = Path(__file__).parent / "data" / "providers" / "argannot" / "ARG-ANNOT_NT_V6_July2019.txt"

# 60 bp synthetic nucleotide sequences, distinct per record (B1 lesson:
# generated from repeat units, never hand-typed).
SEQ_A = "ACGT" * 15
UNIT_B = "ACGTTGCAGG"
SEQ_B = UNIT_B * 6
SEQ_C = "TTAAACCGGG" * 6

EXPECTED_NORMAL = Record(
    db="argannot",
    gene="(AGly)demoA-Ie",
    sequence=SEQ_A,
    accession="FAKE0001:101-709",  # upstream x[1]:x[2]; the trailing length field is ignored
    function=(),
    product="(AGly)demoA-Ie",  # no header description -> gene (save_fasta DESC || ID)
    source_id="(AGly)demoA-Ie:FAKE0001:101-709:609",
)

EXPECTED_BACKSLASH = Record(
    db="argannot",
    gene="(Tet)demoB-Tet",
    sequence=SEQ_B,  # backslash-free: the repair happened before parsing
    accession="FAKE0002:201-812",
    function=(),
    product="(Tet)demoB-Tet",
    source_id="(Tet)demoB-Tet:FAKE0002:201-812:612",
)


@pytest.fixture
def workdir(tmp_path: Path) -> Path:
    """A download workdir holding the fixture under the name the provider's
    source URL would land it as."""
    (tmp_path / "ARG-ANNOT_NT_V6_July2019.txt").write_bytes(FIXTURE.read_bytes())
    return tmp_path


def test_provider_metadata_matches_spec() -> None:
    """Given the argannot PROVIDER, When inspected, Then name, description,
    source URLs, and dbtype match the Wave B7 spec."""
    assert PROVIDER.name == "argannot"
    assert PROVIDER.description == "ARG-ANNOT acquired resistance genes"
    assert PROVIDER.source_urls == (
        "https://web.archive.org/web/20200626214628id_/"
        "https://www.mediterranee-infection.com/wp-content/uploads/2019/09/"
        "ARG-ANNOT_NT_V6_July2019.txt",
    )
    assert PROVIDER.dbtype == "nucl"


def test_transform_yields_exact_records(workdir: Path) -> None:
    """Given the fixture's two well-formed records — one whose raw sequence
    lines are corrupted with trailing and mid-line backslashes — When
    transformed, Then exactly those Records survive in file order: gene is
    the first colon field, accession the x[1]:x[2] join, product the
    description or the gene, source_id the full id token, and the sequence
    contains no backslashes (the upstream repair)."""
    assert tuple(transform(workdir)) == (EXPECTED_NORMAL, EXPECTED_BACKSLASH)


def test_transform_skips_records_with_short_ids(workdir: Path) -> None:
    """Given the fixture also holds a record whose id has only 2
    colon-separated fields, When transformed, Then it is skipped entirely
    (upstream perl would splice an undef into the accession; a typed Record
    refuses to carry that)."""
    records = tuple(transform(workdir))
    assert len(records) == 2
    assert all(record.source_id != "(Fos)demoC-FosA:FAKE0003" for record in records)
