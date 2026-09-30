"""Cluster-kind surface tests: schema selectors, Database.kind discovery,
and `gapit db list --json` kind exposure (stage 1 of the cluster feature).
"""

import json
from collections.abc import Iterable
from pathlib import Path

import pytest
from typer.testing import CliRunner

from gapit.cli import app
from gapit.db import discover_databases
from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record, read_manifest

DATA = Path(__file__).parent / "data" / "cluster"

runner = CliRunner()


def test_schema_features_selector() -> None:
    """Given the new features document, When `gapit schema features`, Then
    its JSON Schema prints with the locus table property."""
    result = runner.invoke(app, ["schema", "features"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert "loci" in payload["properties"]


def test_schema_typing_selector() -> None:
    """Given the new typing document, When `gapit schema typing`, Then its
    JSON Schema prints with the rules/cutoff/margin/fallback properties."""
    result = runner.invoke(app, ["schema", "typing"])
    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    for key in ("rules", "cutoff", "ambiguity_margin", "fallback"):
        assert key in payload["properties"]


def _legacy_dir(datadir: Path, name: str) -> None:
    """An abricate-style database directory: sequences, no manifest."""
    db_dir = datadir / name
    db_dir.mkdir(parents=True)
    (db_dir / "sequences").write_text(
        ">olddb~~~gene1~~~~~~\nACGTACGTACGTACGTACGT\n", encoding="utf-8"
    )


def test_discover_databases_defaults_legacy_dirs_to_gene(tmp_path: Path) -> None:
    """Given a legacy datadir (sequences, no manifest), When discovered,
    Then the Database model reports kind gene."""
    _legacy_dir(tmp_path, "legacy")

    (database,) = discover_databases(tmp_path)

    assert database.kind == "gene"


def test_discover_databases_reports_cluster_kind(tmp_path: Path) -> None:
    """Given a datadir with a cluster-built database, When discovered, Then
    the Database model reports kind cluster."""
    datadir = tmp_path / "dd"
    datadir.mkdir()
    result = runner.invoke(
        app, ["db", "build", "cps", str(DATA / "bakta_style.gbk"), "--datadir", str(datadir)]
    )
    assert result.exit_code == 0, result.stderr
    assert read_manifest(datadir / "cps" / "gapit-manifest.json").kind == "cluster"

    (database,) = discover_databases(datadir)

    assert database.kind == "cluster"


def _syn_provider(name: str) -> Provider:
    """A synthetic registry provider whose upstream is never fetched (the
    db-list tests only need REGISTRY shape; `db build` stays local)."""

    def transform(workdir: Path) -> Iterable[Record]:
        return (
            Record(db=name, gene=record.id, sequence=record.sequence)
            for record in iter_fasta(workdir / "upstream.fa")
        )

    return Provider(
        name=name,
        description="synthetic provider for the cluster-kind tests",
        source_urls=("file://unused",),
        dbtype="nucl",
        transform=transform,
    )


def test_db_list_json_exposes_kind(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given a patched provider registry with one cluster-built provider and
    one legacy installed provider, When `db list --json`, Then every entry
    carries kind — cluster from the manifest, gene for the legacy dir."""
    syn, other = "synamr", "synavail"
    monkeypatch.setattr(
        "gapit.db_ops.REGISTRY", {syn: _syn_provider(syn), other: _syn_provider(other)}
    )
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    built = runner.invoke(
        app, ["db", "build", syn, str(DATA / "bakta_style.gbk"), "--datadir", str(datadir)]
    )
    assert built.exit_code == 0, built.stderr
    _legacy_dir(datadir, other)

    result = runner.invoke(app, ["db", "list", "--datadir", str(datadir), "--json"])

    assert result.exit_code == 0, result.stderr
    document = json.loads(result.stdout)
    by_name = {entry["name"]: entry for entry in document["providers"]}
    assert by_name[syn]["kind"] == "cluster"
    assert by_name[other]["kind"] == "gene"
