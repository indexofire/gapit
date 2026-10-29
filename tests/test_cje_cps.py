"""Bundled cje_cps: the C. jejuni Penner capsule multiplex panel.

The panel is 36 amplicons from Poly et al. 2015 (PLoS ONE 10:e0144349,
CC0) plus an lpxA control. These tests build synthetic carriers from the
wheel's own target sequences, screen them through the real CLI, and pin
the penner_cps calls — including the shared-amplicon profiles (HS5/31,
HS58, HS60, HS32 and the HS4 complex/HS16/HS64 split).
"""

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from gapit.bundled import find_bundled
from gapit.cli import app
from gapit.dbcodec import decode_seqid
from gapit.fasta import iter_fasta

runner = CliRunner()


def _targets() -> dict[str, str]:
    """The bundled target gene -> sequence map (gapit/v1 headers decoded)."""
    bundled = find_bundled("cje_cps")
    assert bundled is not None
    return {
        decode_seqid(record.id, default_db="cje_cps").gene: record.sequence
        for record in iter_fasta(bundled.sequences_path)
    }


def _write_sample(path: Path, genes: list[str]) -> None:
    targets = _targets()
    path.write_text("".join(f">{gene}\n{targets[gene]}\n" for gene in genes), encoding="utf-8")


def _call(sample: Path, datadir: Path, tmp_path: Path) -> str:
    table = tmp_path / f"{sample.stem}.tsv"
    screened = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            "cje_cps",
            "--datadir",
            str(datadir),
            "--output",
            str(table),
            "--quiet",
        ],
    )
    assert screened.exit_code == 0, screened.stderr
    result = runner.invoke(
        app,
        ["typing", str(table), "--datadir", str(datadir), "--format", "json", "--quiet"],
    )
    assert result.exit_code == 0, result.stderr
    document = json.loads(result.stdout)
    return document["files"][0]["phenotypes"]["penner_cps"]["phenotype"]


@pytest.mark.parametrize(
    ("genes", "phenotype"),
    [
        (["lpxA", "HS2"], "HS2"),
        (["lpxA", "HS4A"], "HS4 complex (HS4/13/43/50/62/65)"),
        (["lpxA", "HS4B"], "HS4/13/64 (CG8486)"),
        (["lpxA", "HS4A", "HS4B", "HS52"], "HS16"),
        (["lpxA", "HS4A", "HS4B"], "HS64"),
        (["lpxA", "HS5_31", "HS45"], "HS5/31"),
        (["lpxA", "HS45", "HS60"], "HS60"),
        (["lpxA", "HS15", "HS58"], "HS58"),
        (["lpxA", "HS32", "HS45"], "HS32"),
        (["lpxA", "HS45"], "HS45"),
        (["HS2"], "untypeable (atypical CPS profile or non-jejuni)"),
    ],
)
def test_cje_cps_penner_calls(tmp_path: Path, genes: list[str], phenotype: str) -> None:
    """Given a synthetic carrier of the listed bundled targets, When the
    two-stage CLI runs, Then penner_cps returns the expected call — the
    lpxA-free case exercises the control-gene fallback."""
    datadir = tmp_path / "datadir"
    sample = tmp_path / "sample.fa"
    _write_sample(sample, genes)
    assert _call(sample, datadir, tmp_path) == phenotype


def test_bundled_hs45_alleles_floors_and_screen_defaults() -> None:
    """The curated panel carries the 89-94% HS45 alleles, a 90% default
    identity floor (lpxA/HS12 at 80%), and per-database reads defaults that
    keep MAPQ 20 + a 10% lpxA relative-depth gate."""
    bundled = find_bundled("cje_cps")
    assert bundled is not None
    accessions = {
        decode_seqid(record.id, default_db="cje_cps").accession
        for record in iter_fasta(bundled.sequences_path)
    }
    assert {"HS45_2", "HS45_3", "HS45_4", "HS45_5"} <= accessions
    floors = json.loads((bundled.path / "floors.json").read_text(encoding="utf-8"))
    assert floors["default"] == 90.0
    assert floors["genes"] == {"lpxA": 80.0, "HS12": 80.0}
    screen = json.loads((bundled.path / "screen.json").read_text(encoding="utf-8"))
    assert screen["reads"]["min_mapq"] == 20
    assert screen["reads"]["min_relative_depth"] == 0.10
    assert screen["contigs"]["merge_fragments"] is True


def test_hs45_divergent_allele_is_detected_with_default_floors(tmp_path: Path) -> None:
    """Given the HS45_2 allele added to the panel, When the default screen
    runs (no --minid override), Then the allele clears the 90% floor and the
    panel calls HS45 — the reason the allele was curated."""
    bundled = find_bundled("cje_cps")
    assert bundled is not None
    by_accession = {
        decode_seqid(record.id, default_db="cje_cps").accession: record.sequence
        for record in iter_fasta(bundled.sequences_path)
    }
    sample = tmp_path / "hs45_2.fa"
    sample.write_text(
        f">lpxA\n{by_accession['AL111168.1:251030-251821']}\n>HS45_2\n{by_accession['HS45_2']}\n",
        encoding="utf-8",
    )
    assert _call(sample, tmp_path / "datadir", tmp_path) == "HS45"
