"""CLI integration tests: `gapit db fetch NAME` (provider path) + `db list`.

Offline by construction: ``gapit.cmd_db.REGISTRY`` is monkeypatched down to
two synthetic Providers whose source_urls point at a ``file://`` fasta fixture
in tmp_path, so the full download -> transform -> build pipeline runs against
the pixi env's real makeblastdb/minimap2 without touching the network (the
test_providers_common.py pattern, driven through the CLI surface instead).
"""

import json
from collections.abc import Iterable
from pathlib import Path

import pytest
from pydantic import TypeAdapter
from typer.testing import CliRunner, Result

from gapit.cli import app
from gapit.errors import ErrorEnvelope
from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record

SYN = "synamr"
OTHER = "synavail"
SYN_DESCRIPTION = "synthetic provider for the db fetch CLI tests"

# 60 bp sequences (Wave A3 pinned minimap2 indexing 60-70 bp synthetic seqs).
SEQ_A = "ACGTAG" * 10
SEQ_B = "ACGGCTAGAT" * 6
UPSTREAM_FASTA = f">syn_a first synthetic gene\n{SEQ_A}\n>syn_b second synthetic gene\n{SEQ_B}\n"

runner = CliRunner()
envelope_adapter = TypeAdapter(ErrorEnvelope)


def last_envelope(stderr: str) -> ErrorEnvelope:
    """Parse the last non-empty stderr line as the error envelope."""
    lines = [line for line in stderr.splitlines() if line.strip()]
    assert lines, "expected an error envelope on stderr"
    return envelope_adapter.validate_json(lines[-1])


def syn_transform(workdir: Path) -> Iterable[Record]:
    """Load the downloaded upstream fasta and yield Records (db = SYN)."""
    for fasta in iter_fasta(workdir / "upstream.fa"):
        yield Record(db=SYN, gene=fasta.id, sequence=fasta.sequence)


def patch_registry(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Point ``gapit.cmd_db.REGISTRY`` at two synthetic file:// providers
    (SYN fetchable, OTHER never fetched) and create an empty datadir under
    tmp_path; returns that datadir."""
    source = tmp_path / "upstream.fa"
    source.write_text(UPSTREAM_FASTA, encoding="utf-8")
    url = source.as_uri()
    monkeypatch.setattr(
        "gapit.cmd_db.REGISTRY",
        {
            SYN: Provider(
                name=SYN,
                description=SYN_DESCRIPTION,
                source_urls=(url,),
                dbtype="nucl",
                transform=syn_transform,
            ),
            OTHER: Provider(
                name=OTHER,
                description="available but never fetched in these tests",
                source_urls=(url,),
                dbtype="nucl",
                transform=syn_transform,
            ),
        },
    )
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    return datadir


def fetch(*extra: str) -> Result:
    return runner.invoke(app, ["db", "fetch", *extra])


def test_db_fetch_builds_provider_database(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given a registry provider backed by a file:// fasta, When
    `db fetch NAME --datadir D`, Then exit 0, stdout is a one-line JSON
    receipt (db/records/dbtype/destination) naming D, and the database
    artifacts exist under D (the --datadir flag is respected)."""
    datadir = patch_registry(tmp_path, monkeypatch)
    result = fetch(SYN, "--datadir", str(datadir))
    assert result.exit_code == 0
    db_dir = datadir / SYN
    assert json.loads(result.stdout) == {
        "db": SYN,
        "records": 2,
        "dbtype": "nucl",
        "destination": str(db_dir),
    }
    for name in ("sequences", "sequences.nin", "gapit-manifest.json"):
        assert (db_dir / name).is_file(), name


def test_db_fetch_unknown_provider_is_usage_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a patched registry, When `db fetch` names a provider that is not
    in it, Then exit 2 with a USAGE_ERROR envelope whose message lists the
    available provider names."""
    patch_registry(tmp_path, monkeypatch)
    result = fetch("ncbi")
    assert result.exit_code == 2
    assert result.stdout == ""
    envelope = last_envelope(result.stderr)
    assert envelope.code == "USAGE_ERROR"
    assert "ncbi" in envelope.message
    assert SYN in envelope.message
    assert OTHER in envelope.message


def test_db_fetch_refuses_overwrite_without_force(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given an already-fetched database, When `db fetch NAME` runs again
    without --force, Then exit 4 with a DB_ALREADY_EXISTS envelope and nothing
    on stdout."""
    datadir = patch_registry(tmp_path, monkeypatch)
    assert fetch(SYN, "--datadir", str(datadir)).exit_code == 0
    result = fetch(SYN, "--datadir", str(datadir))
    assert result.exit_code == 4
    assert result.stdout == ""
    envelope = last_envelope(result.stderr)
    assert envelope.code == "DB_ALREADY_EXISTS"
    assert envelope.context["db"] == SYN


def test_db_fetch_force_rebuilds_existing_database(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given an already-fetched database, When `db fetch NAME --force` runs,
    Then exit 0 with a fresh receipt for the same record count."""
    datadir = patch_registry(tmp_path, monkeypatch)
    assert fetch(SYN, "--datadir", str(datadir)).exit_code == 0
    result = fetch(SYN, "--datadir", str(datadir), "--force")
    assert result.exit_code == 0
    assert json.loads(result.stdout)["records"] == 2


def test_db_fetch_creates_missing_nested_datadir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a --datadir whose nested path does not exist yet (fresh
    machine), When `db fetch NAME`, Then exit 0 with the receipt naming the
    created datadir and the database artifacts land under it — not a
    DATADIR_NOT_FOUND refusal (Wave E bootstrap gap)."""
    datadir = patch_registry(tmp_path, monkeypatch)
    nested = datadir / "fresh" / "machine" / "db"
    result = fetch(SYN, "--datadir", str(nested))
    assert result.exit_code == 0
    assert json.loads(result.stdout)["destination"] == str(nested / SYN)
    assert (nested / SYN / "gapit-manifest.json").is_file()


def test_db_list_text_shows_installed_state(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a datadir holding one fetched provider, When `db list`, Then the
    text table shows that provider as installed with its record count and the
    unfetched one as available."""
    datadir = patch_registry(tmp_path, monkeypatch)
    assert fetch(SYN, "--datadir", str(datadir)).exit_code == 0
    result = runner.invoke(app, ["db", "list", "--datadir", str(datadir)])
    assert result.exit_code == 0
    lines = result.stdout.splitlines()
    assert lines[0] == "PROVIDER\tSTATUS\tDBTYPE\tDESCRIPTION"
    assert f"{SYN}\tinstalled (2)\tnucl\t{SYN_DESCRIPTION}" in lines
    assert f"{OTHER}\tavailable\tnucl\tavailable but never fetched in these tests" in lines


def test_db_list_json_emits_gapit_dblist_document(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a datadir holding one fetched provider, When `db list --json`,
    Then stdout parses as a gapit.dblist/1 document whose installed entry
    carries records and whose available entry omits the field."""
    datadir = patch_registry(tmp_path, monkeypatch)
    assert fetch(SYN, "--datadir", str(datadir)).exit_code == 0
    result = runner.invoke(app, ["db", "list", "--datadir", str(datadir), "--json"])
    assert result.exit_code == 0
    document = json.loads(result.stdout)
    assert document["schema"] == "gapit.dblist/1"
    by_name = {entry["name"]: entry for entry in document["providers"]}
    assert by_name[SYN]["installed"] is True
    assert by_name[SYN]["records"] == 2
    assert by_name[SYN]["dbtype"] == "nucl"
    assert by_name[OTHER]["installed"] is False
    assert "records" not in by_name[OTHER]
