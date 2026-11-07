"""Smoke tests for the bundled Vibrio marker/cluster databases."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from gspec.bundled import find_bundled
from gspec.cli import app
from gspec.dbcodec import decode_seqid
from gspec.fasta import iter_fasta

runner = CliRunner()


def _genes(name: str) -> set[str]:
    bundled = find_bundled(name)
    assert bundled is not None
    return {decode_seqid(record.id, name).gene for record in iter_fasta(bundled.sequences_path)}


def _sample(name: str, genes: set[str], destination: Path) -> Path:
    bundled = find_bundled(name)
    assert bundled is not None
    with destination.open("w", encoding="utf-8") as handle:
        for record in iter_fasta(bundled.sequences_path):
            if decode_seqid(record.id, name).gene in genes:
                handle.write(f">{record.id} sampled\n{record.sequence}\n")
    return destination


def test_vpa_pai_has_ten_island_loci() -> None:
    """Given the bundled VPaI cluster database, When inspected, Then it has
    ten VPaI-1..VPaI-10 loci with feature tables and no typing scheme."""
    bundled = find_bundled("vpa_pai")
    assert bundled is not None
    assert bundled.metadata.kind == "cluster"
    assert bundled.typing_path is None
    assert bundled.features_path is not None
    loci = [record.id for record in iter_fasta(bundled.sequences_path)]
    assert loci == [f"VPaI-{number}" for number in range(1, 11)]
    document = bundled.features_path.read_text(encoding="utf-8")
    assert '"id": "VPaI-1"' in document
    assert '"id": "VPaI-10"' in document


def test_vpa_pai_screens_an_island(tmp_path: Path) -> None:
    """Given a bundled VPaI-1 sequence, When screened, Then the cluster engine
    calls VPaI-1 at full coverage."""
    sample = _sample("vpa_pai", {"VPaI-1"}, tmp_path / "vpai1.fna")
    result = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            "vpa_pai",
            "--datadir",
            str(tmp_path / "datadir"),
            "--quiet",
            "--nopath",
        ],
    )
    assert result.exit_code == 0, result.stderr
    row = result.stdout.splitlines()[1].split("\t")
    assert row[1] == "VPaI-1"
    assert row[3] == "100.00"


def test_vch_vf_has_ctx_and_tcpA_markers() -> None:
    """Given the bundled CTX/tcpA gene panel, When inspected, Then the CTX
    core, RS2, ctxB alleles, and the three tcpA markers are all present."""
    genes = _genes("vch_vf")
    assert {
        "ctxA",
        "ctxB_classical",
        "ctxB_eltor",
        "zot",
        "ace",
        "cep",
        "orfU",
        "rstR",
        "rstA",
        "rstB",
        "tcpA_classical",
        "tcpA_eltor",
        "tcpA_o139",
    } == genes


def test_vch_vf_screens_ctx_and_tcpA(tmp_path: Path) -> None:
    """Given CTX/tcpA markers from the bundle, When screened, Then the gene
    engine reports the expected markers at full coverage."""
    wanted = {"ctxA", "ctxB_eltor", "rstR", "tcpA_eltor"}
    sample = _sample("vch_vf", wanted, tmp_path / "vch.fna")
    result = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            "vch_vf",
            "--datadir",
            str(tmp_path / "datadir"),
            "--quiet",
            "--nopath",
        ],
    )
    assert result.exit_code == 0, result.stderr
    rows = [line.split("\t") for line in result.stdout.splitlines()[1:]]
    assert {row[5] for row in rows} == wanted
    assert all(row[9] == "100.00" for row in rows)
