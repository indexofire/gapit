"""CARD provider transform tests (Wave B2) — offline, tarball built at test time.

The committed fixture is plain ``card.json``; each test wraps it into the exact
shape the real download produces: a tar.bz2 (with ``card.json`` at its root)
saved as ``<workdir>/data``, because the CARD source URL
``https://card.mcmaster.ca/latest/data`` has no extension but serves a tarball.
"""

import json
import tarfile
from pathlib import Path

import pytest

from gapit.errors import DatabaseError
from gapit.providers.card import (
    DBTYPE,
    DESCRIPTION,
    NAME,
    PROVIDER,
    SOURCE_URLS,
    transform,
)
from gapit.records import Record

_FIXTURE = Path(__file__).parent / "data" / "providers" / "card" / "card.json"


def _workdir(card_payload: dict[str, object], tmp_path: Path, arcname: str = "card.json") -> Path:
    """Given a card dict, wrap it into <workdir>/data (member named ``arcname``)."""
    card_path = tmp_path / "card.json"
    card_path.write_text(json.dumps(card_payload), encoding="utf-8")
    workdir = tmp_path / "work"
    workdir.mkdir()
    with tarfile.open(workdir / "data", "w:bz2") as tar:
        tar.add(card_path, arcname=arcname)
    return workdir


def _fixture_payload() -> dict[str, object]:
    """The committed fixture, parsed (fresh copy per call, safe to mutate)."""
    payload: dict[str, object] = json.loads(_FIXTURE.read_text(encoding="utf-8"))
    return payload


def _fixture_records(tmp_path: Path) -> tuple[Record, ...]:
    """Run the transform over the committed fixture (through its tarball shape)."""
    return tuple(transform(_workdir(_fixture_payload(), tmp_path)))


def test_transform_when_typical_model_yields_exact_record(tmp_path: Path) -> None:
    # Given: the committed CARD fixture, wrapped as the real tarball download
    # When: the card transform runs
    records = {record.gene: record for record in _fixture_records(tmp_path)}
    # Then: the plain plus-strand model is exact field-for-field — the
    # SORTED-first sequence key wins (not insertion-first; the decoy
    # homolog-seq-2 sorts after homolog-seq-1), and the lowercase 'g' passes
    # through verbatim (normalization is fetch_provider's job, not the
    # transform's).
    oxa1 = records["OXA-1"]
    assert oxa1.db == "card"
    assert oxa1.gene == "OXA-1"
    assert oxa1.sequence == "ATGgCAAAAACTCGACTACGCACTTCGCTAAGCTTGGCTTACGG"
    assert oxa1.accession == "NC_000913.3:100-900"
    assert oxa1.function == ("penam", "cephalosporin")
    assert oxa1.product == "class A beta-lactamase OXA-1"
    assert oxa1.source_id == "1000001"


def test_transform_when_minus_strand_swaps_coordinates(tmp_path: Path) -> None:
    # Given: the fixture's strand '-' model (model_name "Qnr  B19", two spaces)
    # When: the card transform runs
    records = {record.gene: record for record in _fixture_records(tmp_path)}
    # Then: coordinates swap to accession:fmax-fmin, and the model_name
    # whitespace RUN collapses to a single underscore (perl s/\s+/_/g).
    qnr = records["Qnr_B19"]
    assert qnr.accession == "NZ_CP044330.1:990-210"
    assert qnr.gene == "Qnr_B19"
    assert qnr.product == "quinolone resistance protein"
    assert qnr.function == ()


def test_transform_when_multiple_drug_classes(tmp_path: Path) -> None:
    # Given: a model with three Drug Class categories (one without the
    # " antibiotic" suffix, one with internal spaces) and no ARO_description
    # When: the card transform runs
    records = {record.gene: record for record in _fixture_records(tmp_path)}
    # Then: " antibiotic" is stripped (once), whitespace becomes '_' per
    # character, non-Drug-Class categories are filtered out, insertion order
    # is preserved (sorting is fetch_provider's job), and the missing
    # ARO_description falls back to ARO_accession (perl ||).
    aac = records["AAC3-IV"]
    assert aac.function == ("macrolide", "quaternary_ammonium_compound", "phenicol")
    assert aac.product == "3002870"
    assert aac.source_id == "3002870"
    assert aac.accession == "M81112.1:5-805"


def test_transform_when_non_homolog_models_are_skipped(tmp_path: Path) -> None:
    # Given: the fixture contains a "protein variant model" (which even carries
    # model_param.snp — the model_type filter must run BEFORE the snp error)
    # plus the "_version" string entry
    # When: the card transform runs
    records = _fixture_records(tmp_path)
    # Then: exactly the three protein homolog models remain, in file order.
    assert [record.gene for record in records] == ["OXA-1", "Qnr_B19", "AAC3-IV"]


def test_transform_when_snp_homolog_model_raises(tmp_path: Path) -> None:
    # Given: a protein homolog model carrying model_param.snp
    payload = _fixture_payload()
    payload["1000005"] = {
        "model_name": "mcr-9 homolog",
        "model_type": "protein homolog model",
        "ARO_accession": "1000005",
        "model_param": {"snp": {"original": "A", "position": 21, "change": "G"}},
        "model_sequences": {"sequence": {}},
    }
    # When: the transform is consumed
    # Then: upstream's err() surfaces as a typed DatabaseError with the RAW
    # model_name in context (before any whitespace substitution).
    with pytest.raises(DatabaseError) as raised:
        list(transform(_workdir(payload, tmp_path)))
    assert raised.value.code == "PROVIDER_INVALID"
    assert raised.value.context == {"model": "mcr-9 homolog"}


def test_transform_when_members_are_dot_slash_prefixed(tmp_path: Path) -> None:
    # Given: the re-tarred card.mcmaster.ca archive prefixes members with
    # './' (Wave E: the exact-name extractfile("card.json") lookup missed)
    # When: the transform runs over a './card.json' tarball
    records = tuple(transform(_workdir(_fixture_payload(), tmp_path, arcname="./card.json")))
    # Then: the member is matched by root-normalized name and the exact
    # records still come out, in file order.
    assert [record.gene for record in records] == ["OXA-1", "Qnr_B19", "AAC3-IV"]


def test_transform_when_card_json_member_is_missing_raises(tmp_path: Path) -> None:
    # Given: a tarball whose only member is not card.json
    other = tmp_path / "other.json"
    other.write_text(json.dumps(_fixture_payload()), encoding="utf-8")
    workdir = tmp_path / "work"
    workdir.mkdir()
    with tarfile.open(workdir / "data", "w:bz2") as tar:
        tar.add(other, arcname="./other.json")
    # When: the transform is consumed
    # Then: the missing member surfaces as a typed PROVIDER_INVALID (exit 4)
    # naming the archive and the expected member — not an untyped KeyError
    # escaping as UNEXPECTED/exit 1 (Wave E error-contract violation).
    with pytest.raises(DatabaseError) as raised:
        list(transform(workdir))
    error = raised.value
    assert error.code == "PROVIDER_INVALID"
    assert error.exit_code == 4
    assert error.context == {"archive": str(workdir / "data"), "expected": "card.json"}


def test_provider_metadata_matches_the_b0_contract() -> None:
    # Given/When: the card provider module's exported values
    # Then: they match the pinned provider metadata and Provider wiring.
    assert NAME == "card"
    assert DESCRIPTION == "CARD protein homolog resistance models"
    assert SOURCE_URLS == ("https://card.mcmaster.ca/latest/data",)
    assert DBTYPE == "nucl"
    assert PROVIDER.name == NAME
    assert PROVIDER.description == DESCRIPTION
    assert PROVIDER.source_urls == SOURCE_URLS
    assert PROVIDER.dbtype == DBTYPE
    assert PROVIDER.transform is transform
