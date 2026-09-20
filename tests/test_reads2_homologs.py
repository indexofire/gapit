"""Family-splitting regression tests (the gapit.reads/2 money tests).

Fixture (tests/data/reads2_db/homologs + tests/data/reads2; seeded recipe in
the generator's docstring, probe-verified on minimap2 2.31): the db holds a
600 nt geneA (truly carried) plus a 400 nt partial homolog geneB (~85%
identical, absent); the sample also carries a novel B-like allele C (~95% to
geneB) whose reads pile onto geneB as ~90%-identity primary alignments.
Unfiltered reads mode therefore calls BOTH db genes present — the documented
family-splitting/over-call failure (KP benchmark: sr 11 true names vs ONT 140
vs blastn 18). ``--min-identity 95`` drops the ~90% alignments on geneB and
keeps the ~100% (sr) / ~98% (ONT) ones on geneA: only the true gene survives.

The sr leg calibrates --min-breadth 80: sr soft-clipping and the no-tag PAF
geometry cap geneB at ~89.8% breadth (docs/reads.md documents calibrating
breadth for short targets); the ONT leg over-calls at the default 90.
"""

import shutil
from pathlib import Path

import pytest
from pydantic import TypeAdapter
from typer.testing import CliRunner

from gapit.cli import app
from gapit.formats.json import Reads2Document, ReadsDocument

READS2_DB = Path(__file__).parent / "data" / "reads2_db"
READS2 = Path(__file__).parent / "data" / "reads2"

runner = CliRunner()
reads1_adapter = TypeAdapter(ReadsDocument)
reads2_adapter = TypeAdapter(Reads2Document)


@pytest.fixture()
def datadir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    target = tmp_path / "datadir"
    shutil.copytree(READS2_DB, target)
    monkeypatch.chdir(READS2)
    return target


def _screen(datadir: Path, fq: str, read_type: str, extra: list[str]) -> str:
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            fq,
            "--db",
            "homologs",
            "--datadir",
            str(datadir),
            "--read-type",
            read_type,
            "--quiet",
            *extra,
        ],
    )
    assert result.exit_code == 0, result.stderr
    return result.stdout


def screen1(datadir: Path, fq: str, read_type: str, *extra: str) -> ReadsDocument:
    return reads1_adapter.validate_json(_screen(datadir, fq, read_type, list(extra)))


def screen2(datadir: Path, fq: str, read_type: str, *extra: str) -> Reads2Document:
    return reads2_adapter.validate_json(
        _screen(datadir, fq, read_type, [*extra, "--min-identity", "95"])
    )


def present_genes(document: ReadsDocument | Reads2Document) -> set[str]:
    (first,) = document.files
    return {entry.gene for entry in first.genes if entry.present}


def gene_names(document: ReadsDocument | Reads2Document) -> set[str]:
    (first,) = document.files
    return {entry.gene for entry in first.genes}


def test_sr_homolog_split_unfiltered_calls_both(datadir: Path) -> None:
    """Given short-read fixtures over a two-homolog db, When screened
    unfiltered at a calibrated breadth floor, Then BOTH genes are called
    present (the reproduced family-splitting over-call)."""
    document = screen1(datadir, "sr_homologs.fq", "sr", "--min-breadth", "80")
    assert document.schema_name == "gapit.reads/1"
    assert present_genes(document) == {"geneA", "geneB"}


def test_sr_min_identity_95_calls_only_true_gene(datadir: Path) -> None:
    """Given the same fixtures with --min-identity 95, When screened, Then
    only geneA is present and geneB vanishes from the gene list (its ~90%
    alignments are all dropped; zero-read genes are omitted)."""
    document = screen2(datadir, "sr_homologs.fq", "sr", "--min-breadth", "80")
    assert document.schema_name == "gapit.reads/2"
    assert gene_names(document) == {"geneA"}
    assert present_genes(document) == {"geneA"}


def test_sr_true_gene_keeps_full_identity(datadir: Path) -> None:
    """Given the filtered sr run, When inspected, Then geneA's
    mean_identity_pct is ~100 (clean reads) and its breadth still passes the
    presence threshold."""
    document = screen2(datadir, "sr_homologs.fq", "sr", "--min-breadth", "80")
    (first,) = document.files
    (entry,) = first.genes
    assert entry.gene == "geneA"
    assert entry.mean_identity_pct == 100.0
    assert entry.breadth_pct >= 90.0
    assert entry.present is True


def test_ont_homolog_split_unfiltered_calls_both(datadir: Path) -> None:
    """Given ONT-style noisy fixtures, When screened unfiltered with
    map-ont at the default breadth, Then BOTH genes are called present (the
    ONT over-call family from the KP benchmark)."""
    document = screen1(datadir, "ont_homologs.fq", "map-ont")
    assert present_genes(document) == {"geneA", "geneB"}


def test_ont_min_identity_95_calls_only_true_gene(datadir: Path) -> None:
    """Given the ONT fixtures with --min-identity 95, When screened, Then
    only geneA survives (its noisy alignments stay ~98%, geneB's ~89%
    alignments are all dropped and the gene is omitted)."""
    document = screen2(datadir, "ont_homologs.fq", "map-ont")
    (first,) = document.files
    assert gene_names(document) == {"geneA"}
    (entry,) = first.genes
    assert entry.mean_identity_pct >= 95.0
    assert entry.breadth_pct == 100.0
    assert entry.present is True
