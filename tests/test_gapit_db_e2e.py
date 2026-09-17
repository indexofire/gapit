"""Offline end-to-end: gapit/v1 native databases through the real screening
path.

Builds a tiny database from synthetic records via dbbuild (records.jsonl ->
sequences + BLAST index + .mmi + manifest, written last), then screens a
synthetic contig (contig mode, --format json) and synthetic 100 bp reads
(reads mode, default json — the run goes through the persisted .mmi, whose
manifest certifies the same minimap2 as the installed one). Uses the pixi
env's real makeblastdb/blastn/minimap2, same pattern as
test_screen_integration.py and test_dbbuild.py; everything is synthetic and
local to tmp_path, no network.
"""

import random
from pathlib import Path

import pytest
from pydantic import TypeAdapter
from typer.testing import CliRunner

from gapit.cli import app
from gapit.dbbuild import build_database
from gapit.formats.json import ReadsDocument, ReportDocument
from gapit.records import Record, write_records

DB = "tinygapit"
SEQ_LEN = 240


def _synthetic_dna(length: int, seed: int) -> str:
    """Deterministic pseudo-random ACGT (random.Random(int) is reproducible
    across runs/platforms). Non-repetitive on purpose: minimap2 -x sr clips
    alignments hard on repetitive targets (repeat-unit fixtures capped
    breadth at ~71%); distinct seeds share no 11-mer in practice, so the two
    genes cannot cross-match under blastn or minimap2 either."""
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


SEQ_DEMO = _synthetic_dna(SEQ_LEN, seed=1)
SEQ_OTHER = _synthetic_dna(SEQ_LEN, seed=2)

RECORDS: tuple[Record, ...] = (
    Record(
        db=DB,
        gene="de|mo%A=1",
        sequence=SEQ_DEMO,
        accession="SYN-DEMO-001",
        function=("ampicillin", "gentamicin"),
        product="demo beta-lactamase, variant A",
    ),
    Record(
        db=DB,
        gene="otherv2",
        sequence=SEQ_OTHER,
        accession="SYN-DEMO-002",
        function=("tetracycline",),
        product="demo efflux pump",
    ),
)

runner = CliRunner()
report_adapter = TypeAdapter(ReportDocument)
reads_adapter = TypeAdapter(ReadsDocument)


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    """A fully built gapit/v1 datadir: all four artifacts, manifest last."""
    db_dir = tmp_path / "datadir" / DB
    db_dir.mkdir(parents=True)
    write_records(RECORDS, db_dir / "records.jsonl")
    build_database(
        db_dir, name=DB, dbtype="nucl", source_urls=(), fetched_at="2026-09-17T00:00:00Z"
    )
    return tmp_path / "datadir"


def _tiled_reads(tmp_path: Path) -> Path:
    """Three overlapping 100 bp reads exactly tiling the demo gene (the sr
    preset soft-clips a few nt per alignment end; 30 bp overlaps keep the
    union breadth safely above the 90% presence threshold)."""
    spans = ((0, 100), (70, 170), (140, SEQ_LEN))
    lines = "".join(
        f"@read{i}\n{SEQ_DEMO[start:stop]}\n+\n{'I' * (stop - start)}\n"
        for i, (start, stop) in enumerate(spans)
    )
    reads = tmp_path / "demo.fq"
    reads.write_text(lines, encoding="utf-8")
    return reads


def test_screen_json_decodes_native_headers(datadir: Path, tmp_path: Path) -> None:
    """Given the built native db and a contig that IS the metachar gene, When
    screened with --format json, Then exactly one hit decodes gene/accession/
    resistance ('a;b' display form), the product lost its id prefix, and the
    second gene is absent."""
    query = tmp_path / "contig.fa"
    query.write_text(f">contig1\n{SEQ_DEMO}\n", encoding="utf-8")
    result = runner.invoke(
        app,
        ["screen", str(query), "--db", DB, "--datadir", str(datadir), "--format", "json"],
    )
    assert result.exit_code == 0
    document = report_adapter.validate_json(result.stdout)
    assert document.schema_name == "gapit.report/1"
    assert document.params.db == DB
    (entry,) = document.files
    assert entry.file == str(query)
    (hit,) = entry.hits
    assert hit.gene == "de|mo%A=1"
    assert hit.database == DB
    assert hit.accession == "SYN-DEMO-001"
    assert hit.resistance == "ampicillin;gentamicin"
    assert hit.product == "demo beta-lactamase variant A"
    assert hit.identity_pct == 100.0
    assert hit.coverage_pct == 100.0


def test_reads_screen_uses_persisted_mmi(datadir: Path, tmp_path: Path) -> None:
    """Given the built db — its sequences.mmi is certified by a manifest that
    matches the installed minimap2, so reads mode screens against the .mmi —
    and reads tiling the demo gene, When reads-screened, Then exactly that
    gene decodes and is present."""
    assert (datadir / DB / "sequences.mmi").is_file()
    result = runner.invoke(
        app,
        ["screen", "--r1", str(_tiled_reads(tmp_path)), "--db", DB, "--datadir", str(datadir)],
    )
    assert result.exit_code == 0
    document = reads_adapter.validate_json(result.stdout)
    assert document.schema_name == "gapit.reads/1"
    (entry,) = document.files[0].genes
    assert entry.gene == "de|mo%A=1"
    assert entry.database == DB
    assert entry.accession == "SYN-DEMO-001"
    assert entry.resistance == "ampicillin;gentamicin"
    assert entry.product == "demo beta-lactamase, variant A"
    assert entry.present is True
    assert entry.reads_mapped == 3


def test_reads_screen_falls_back_to_fasta_without_manifest(datadir: Path, tmp_path: Path) -> None:
    """Given the same db with its manifest deleted (mmi no longer verifiable),
    When reads-screened, Then screening still succeeds against the FASTA —
    a silent performance fallback, not an error — with the same decoded call."""
    (datadir / DB / "gapit-manifest.json").unlink()
    result = runner.invoke(
        app,
        ["screen", "--r1", str(_tiled_reads(tmp_path)), "--db", DB, "--datadir", str(datadir)],
    )
    assert result.exit_code == 0
    document = reads_adapter.validate_json(result.stdout)
    (entry,) = document.files[0].genes
    assert entry.gene == "de|mo%A=1"
    assert entry.present is True
