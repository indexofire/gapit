"""Kaptive cluster-provider suites (stage 3): `gapit db fetch` for the seven
official-keyword providers.

Offline by construction: the REGISTRY is monkeypatched so every kaptive
provider's source URL points at a file:// copy of the kaptive-style GBK
fixture — the full download -> parse -> cluster-build pipeline runs against
the pixi env's real makeblastdb/minimap2 with zero network (the
test_cli_db_fetch.py pattern). The pinned upstream URLs are checked
against the real REGISTRY in the first test (no requests are made).
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

# The seven official Kaptive v3 install keywords, each sourcing the raw
# GenBank file from main of its actively curated per-species repo.
EXPECTED_SOURCES: dict[str, tuple[str, str]] = {
    "kpsc_k": (
        "https://raw.githubusercontent.com/klebgenomics/"
        "KpSC_surface_antigen_loci/main/"
        "Klebsiella_pneumoniae_Species_Complex_K.gbk",
        "klebgenomics/KpSC_surface_antigen_loci",
    ),
    "kpsc_o": (
        "https://raw.githubusercontent.com/klebgenomics/"
        "KpSC_surface_antigen_loci/main/"
        "Klebsiella_pneumoniae_Species_Complex_O.gbk",
        "klebgenomics/KpSC_surface_antigen_loci",
    ),
    "kosc_k": (
        "https://raw.githubusercontent.com/klebgenomics/"
        "KoSC-surface-antigen-loci/main/"
        "Klebsiella_oxytoca_Species_Complex_K_locus_database.gbk",
        "klebgenomics/KoSC-surface-antigen-loci",
    ),
    "kosc_o": (
        "https://raw.githubusercontent.com/klebgenomics/"
        "KoSC-surface-antigen-loci/main/"
        "Klebsiella_oxytoca_Species_Complex_O_locus_database.gbk",
        "klebgenomics/KoSC-surface-antigen-loci",
    ),
    "ab_k": (
        "https://raw.githubusercontent.com/johannajkenyon/"
        "Abaumannii_surface_polysaccharide_loci/main/Acinetobacter_baumannii_K.gbk",
        "johannajkenyon/Abaumannii_surface_polysaccharide_loci",
    ),
    "ab_o": (
        "https://raw.githubusercontent.com/johannajkenyon/"
        "Abaumannii_surface_polysaccharide_loci/main/Acinetobacter_baumannii_OC.gbk",
        "johannajkenyon/Abaumannii_surface_polysaccharide_loci",
    ),
    "ecoli_kps": (
        "https://raw.githubusercontent.com/rgladstone/EC-K-typing/main/"
        "EC-K-typing_group2and3.gbk",
        "rgladstone/EC-K-typing",
    ),
}
KAPTIVE_PROVIDERS = tuple(EXPECTED_SOURCES)


@pytest.fixture()
def kaptive_registry(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Point the seven kaptive providers at the file:// fixture; returns an
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


def test_real_registry_carries_the_seven_kaptive_providers() -> None:
    """Given the shipped provider registry, When inspected, Then the seven
    official-keyword kaptive names are the ONLY cluster providers, each
    with its verified raw main URL, upstream-repo + Wyres provenance, and
    no snapshot (nothing bundled)."""
    assert {n for n, p in REGISTRY.items() if p.kind == "cluster"} == set(KAPTIVE_PROVIDERS)
    for name, (url, repo) in EXPECTED_SOURCES.items():
        provider = REGISTRY[name]
        assert isinstance(provider, ClusterProvider), name
        assert provider.kind == "cluster"
        assert provider.dbtype == "nucl"
        assert provider.snapshot is None
        assert provider.source_urls == (url,), name
        assert provider.license == "GPL-3.0 (database content)"
        assert repo in provider.note, name
        assert "Wyres" in provider.note, name


def test_fetch_builds_cluster_database(
    kaptive_registry: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the file:// kaptive fixture, When `db fetch kpsc_k`, Then
    exit 0 with a receipt counting the fixture's loci and the database
    carries every cluster artifact — locus FASTA, features.json, index,
    manifest kind cluster — but NO records.jsonl and NO typing.json
    (phenotype stays null until a typing spec exists)."""
    result = runner.invoke(app, ["db", "fetch", "kpsc_k", "--datadir", str(kaptive_registry)])
    assert result.exit_code == 0, result.stderr
    db_dir = kaptive_registry / "kpsc_k"
    assert json.loads(result.stdout) == {
        "db": "kpsc_k",
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
    result = runner.invoke(app, ["db", "fetch", "kpsc_o", "--datadir", str(kaptive_registry)])
    assert result.exit_code == 0, result.stderr
    manifest = read_manifest(kaptive_registry / "kpsc_o" / "gapit-manifest.json")
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
            "ecoli_kps": ClusterProvider(
                name="ecoli_kps",
                description="screening fixture stand-in",
                source_urls=(url,),
                license="GPL-3.0 (database content)",
                note="fixture note",
            )
        },
    )
    assert (
        runner.invoke(app, ["db", "fetch", "ecoli_kps", "--datadir", str(datadir)]).exit_code == 0
    )
    locus = "".join(
        line.strip()
        for line in (datadir / "ecoli_kps" / "sequences").read_text(encoding="utf-8").splitlines()
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
            "ecoli_kps",
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
    args = ["db", "fetch", "ab_k", "--datadir", str(kaptive_registry)]
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
            app, ["db", "fetch", "kosc_k", "--datadir", str(kaptive_registry)]
        ).exit_code
        == 0
    )
    result = runner.invoke(app, ["db", "list", "--datadir", str(kaptive_registry), "--json"])
    assert result.exit_code == 0, result.stderr
    entries = {entry["name"]: entry for entry in json.loads(result.stdout)["providers"]}
    assert entries["kosc_k"]["kind"] == "cluster"
    assert entries["kosc_k"]["installed"] is True
    assert entries["kosc_k"]["records"] == 2
    assert entries["kosc_o"]["kind"] == "cluster"
    assert entries["kosc_o"]["installed"] is False
