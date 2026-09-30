"""Kaptive cluster-provider suites (stage 3): `gapit db fetch kaptive_*`.

Offline by construction: the REGISTRY is monkeypatched so every kaptive
provider's source URL points at a file:// copy of the kaptive-style GBK
fixture — the full download -> parse -> cluster-build pipeline runs against
the pixi env's real makeblastdb/minimap2 with zero network (the
test_cli_db_fetch.py pattern). The pinned upstream URLs are separate
one-request integration checks marked network.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from gapit.cli import app
from gapit.providers import REGISTRY
from gapit.providers.cluster_common import ClusterProvider
from gapit.records import read_manifest

DATA = Path(__file__).parent / "data" / "cluster"
FIXTURE = DATA / "kaptive_style.gbk"  # two kaptive-style loci: KL106, KL107
runner = CliRunner()

KAPTIVE_PROVIDERS = ("kaptive_ak", "kaptive_k", "kaptive_o", "kaptive_oc")


@pytest.fixture()
def kaptive_registry(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Point the four kaptive providers at the file:// fixture; returns an
    empty datadir."""
    url = FIXTURE.as_uri()
    monkeypatch.setattr(
        "gapit.db_ops.REGISTRY",
        {
            name: ClusterProvider(
                name=name,
                description=f"fixture {name}",
                source_urls=(url,),
                license="GPL-3.0 (database content)",
                note="fixture note",
            )
            for name in KAPTIVE_PROVIDERS
        },
    )
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    return datadir


def test_real_registry_carries_the_four_kaptive_providers() -> None:
    """Given the shipped provider registry, When inspected, Then the four
    kaptive names are present as cluster providers with the verified v2.0.9
    raw URLs and GPL-3.0 provenance — and no snapshot (nothing bundled)."""
    for name in KAPTIVE_PROVIDERS:
        provider = REGISTRY[name]
        assert isinstance(provider, ClusterProvider), name
        assert provider.kind == "cluster"
        assert provider.dbtype == "nucl"
        assert provider.snapshot is None
        (url,) = provider.source_urls
        assert url.startswith(
            "https://raw.githubusercontent.com/klebgenomics/Kaptive/v2.0.9/reference_database/"
        ), url
        assert url.endswith(".gbk")
        assert provider.license == "GPL-3.0 (database content)"
        assert "Wyres" in provider.note


def test_fetch_builds_cluster_database(
    kaptive_registry: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the file:// kaptive fixture, When `db fetch kaptive_k`, Then
    exit 0 with a receipt counting the fixture's loci and the database
    carries every cluster artifact — locus FASTA, features.json, index,
    manifest kind cluster — but NO records.jsonl and NO typing.json
    (phenotype stays null until a typing spec exists)."""
    result = runner.invoke(app, ["db", "fetch", "kaptive_k", "--datadir", str(kaptive_registry)])
    assert result.exit_code == 0, result.stderr
    db_dir = kaptive_registry / "kaptive_k"
    assert json.loads(result.stdout) == {
        "db": "kaptive_k",
        "records": 2,
        "dbtype": "nucl",
        "destination": str(db_dir),
    }
    headers = [
        line
        for line in (db_dir / "sequences").read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    ]
    assert headers == [">KL106", ">KL107"]
    assert (db_dir / "features.json").is_file()
    assert (db_dir / "sequences.nin").is_file()
    assert not (db_dir / "records.jsonl").exists()
    assert not (db_dir / "typing.json").exists()


def test_fetch_manifest_carries_provenance(kaptive_registry: Path) -> None:
    """Given a fetched kaptive database, When its manifest is read, Then it
    declares kind cluster, the upstream source URL, the database-content
    license, and the citation note."""
    result = runner.invoke(app, ["db", "fetch", "kaptive_o", "--datadir", str(kaptive_registry)])
    assert result.exit_code == 0, result.stderr
    manifest = read_manifest(kaptive_registry / "kaptive_o" / "gapit-manifest.json")
    assert manifest.kind == "cluster"
    assert manifest.header_format == "plain"
    assert manifest.n_records == 2
    assert manifest.source_urls == (FIXTURE.as_uri(),)
    assert manifest.license == "GPL-3.0 (database content)"
    assert manifest.note == "fixture note"


def test_fetched_kaptive_db_screens_as_cluster(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the fetched fixture db and a contig extracted from one locus,
    When screened, Then the untyped cluster engine calls the locus (locus
    call only — phenotype '-', no typing model installed). The source is
    screening.gbk: its 1200 bp loci clear minimap2's asm20 thresholds while
    kaptive_style.gbk's 160 bp toy loci do not."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    url = (DATA / "screening.gbk").as_uri()
    monkeypatch.setattr(
        "gapit.db_ops.REGISTRY",
        {
            "kaptive_oc": ClusterProvider(
                name="kaptive_oc",
                description="screening fixture stand-in",
                source_urls=(url,),
                license="GPL-3.0 (database content)",
                note="fixture note",
            )
        },
    )
    assert (
        runner.invoke(app, ["db", "fetch", "kaptive_oc", "--datadir", str(datadir)]).exit_code == 0
    )
    locus = "".join(
        line.strip()
        for line in (datadir / "kaptive_oc" / "sequences").read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    )
    sample = tmp_path / "sample.fa"
    sample.write_text(f">ctg\n{locus}\n", encoding="utf-8")
    result = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            "kaptive_oc",
            "--datadir",
            str(datadir),
            "--quiet",
        ],
    )
    assert result.exit_code == 0, result.stderr
    header, row = result.stdout.splitlines()
    assert header == "FILE\tBEST_LOCUS\tTYPE\tCOVERAGE\tIDENTITY\tPRESENT\tPARTIAL\tMISSING_IDS"
    assert row.startswith(f"{sample}\tlocusA\tKL101\t100.00\t100.00\t3\t0\t-")


def test_fetch_refuses_overwrite_and_force_rebuilds(kaptive_registry: Path) -> None:
    """Given an already-fetched kaptive database, When fetched again without
    --force, Then DB_ALREADY_EXISTS; with --force the rebuild succeeds."""
    args = ["db", "fetch", "kaptive_ak", "--datadir", str(kaptive_registry)]
    assert runner.invoke(app, args).exit_code == 0
    again = runner.invoke(app, args)
    assert again.exit_code == 4
    assert "DB_ALREADY_EXISTS" in again.stderr
    forced = runner.invoke(app, [*args, "--force"])
    assert forced.exit_code == 0, forced.stderr


def test_db_list_marks_kaptive_providers_cluster_kind(kaptive_registry: Path) -> None:
    """Given the kaptive registry with one database fetched, When
    `db list --json`, Then fetched and unfetched kaptive entries alike carry
    kind cluster (gene providers stay gene)."""
    assert (
        runner.invoke(
            app, ["db", "fetch", "kaptive_k", "--datadir", str(kaptive_registry)]
        ).exit_code
        == 0
    )
    result = runner.invoke(app, ["db", "list", "--datadir", str(kaptive_registry), "--json"])
    assert result.exit_code == 0, result.stderr
    entries = {entry["name"]: entry for entry in json.loads(result.stdout)["providers"]}
    assert entries["kaptive_k"]["kind"] == "cluster"
    assert entries["kaptive_k"]["installed"] is True
    assert entries["kaptive_k"]["records"] == 2
    assert entries["kaptive_o"]["kind"] == "cluster"
    assert entries["kaptive_o"]["installed"] is False
