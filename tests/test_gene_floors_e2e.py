"""E2E: per-gene identity floors through the real CLI (db-driven opt-in).

Leg 1 — the homologs fixture (tests/data/reads2_db + tests/data/reads2): a
gapit.floors/1 sidecar DROPPED INTO the db directory flips the documented
family-splitting over-call exactly like ``--min-identity 95`` would, but
only for the floored gene; removing the sidecar restores the byte-identical
floorless table (nothing about the CLI changed — the database did).

Leg 2 — the motivating SPATE case on the bundled ecoli_dec: a pic panel
sequence mutated to ~87% identity (26ECO0084: 97.59% breadth at 86.6%
identity) is called present without floors and vanishes under the shipped
{pic: 90} floor, while the unfloored true alleles in the same assembly
(uidA, aggR) stay present in both runs.
"""

import random
import shutil
from pathlib import Path

from typer.testing import CliRunner

from gapit.cli import app
from gapit.fasta import iter_fasta

READS2_DB = Path(__file__).parent / "data" / "reads2_db"
READS2 = Path(__file__).parent / "data" / "reads2"
BUNDLE = Path(__file__).parent.parent / "src" / "gapit" / "data" / "dbs" / "ecoli_dec"

runner = CliRunner()


def screen_reads_tsv(datadir: Path, query: Path, db: str, *extra: str) -> str:
    """One reads-mode screen (FASTA input content-detects to map-ont),
    default TSV, quiet; returns stdout for byte comparison."""
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            str(query),
            "--db",
            db,
            "--datadir",
            str(datadir),
            "--quiet",
            *extra,
        ],
    )
    assert result.exit_code == 0, result.stderr
    return result.stdout


def tsv_genes(stdout: str) -> set[str]:
    """GENE (column 2) of every data row."""
    return {line.split("\t")[1] for line in stdout.splitlines() if not line.startswith("#")}


# ------------------------------------------------- leg 1: homologs fixture --


def test_homologs_sidecar_flips_only_the_floored_gene(tmp_path: Path) -> None:
    """Given the two-homolog db and sr fixtures, When the floorless screen
    runs, then a {geneB: 95} floors.json is dropped into the db directory
    and the screen reruns, Then geneB vanishes while geneA stays; removing
    the sidecar restores the byte-identical floorless table."""
    datadir = tmp_path / "datadir"
    shutil.copytree(READS2_DB / "homologs", datadir / "homologs")
    floors = datadir / "homologs" / "floors.json"
    query = READS2 / "sr_homologs.fq"

    before = screen_reads_tsv(datadir, query, "homologs", "--min-breadth", "80")
    assert tsv_genes(before) == {"geneA", "geneB"}  # the reproduced over-call

    floors.write_text(
        '{"schema": "gapit.floors/1", "default": null, "genes": {"geneB": 95.0}}',
        encoding="utf-8",
    )
    floored = screen_reads_tsv(datadir, query, "homologs", "--min-breadth", "80")
    assert tsv_genes(floored) == {"geneA"}

    floors.unlink()
    stripped = screen_reads_tsv(datadir, query, "homologs", "--min-breadth", "80")
    assert stripped == before


# --------------------------------------------------- leg 2: the pic case --


def _mutate(sequence: str, divergence: float, seed: int) -> str:
    """Deterministically substitute ~`divergence` of the bases (each to a
    different base) — a synthetic SPATE homolog."""
    rng = random.Random(seed)
    return "".join(
        base if rng.random() >= divergence else next(b for b in "ACGT" if b != base)
        for base in sequence
    )


def test_pic_homolog_is_absent_with_floors_present_without(tmp_path: Path) -> None:
    """Given the bundled ecoli_dec (ships floors.json {pic: 90}) and an
    assembly of a ~87%-identity pic homolog beside true uidA/aggR alleles,
    When screened in reads mode, Then pic is absent (sub-floor alignments
    dropped, breadth collapsed) while uidA/aggR stay present; after a
    floors-stripped rebuild (sidecar removed) the same assembly calls pic
    present — breadth-only presence was the over-call, the database-side
    floor is the fix."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    # Materialize the bundle by screening any panel-derived sample once.
    true_seqs: list[tuple[str, str]] = []
    pic_id: str | None = None
    pic_seq: str | None = None
    for record in iter_fasta(BUNDLE / "sequences"):
        gene = record.id.split("~~~")[1]
        if gene == "pic" and pic_seq is None:
            pic_id, pic_seq = record.id, record.sequence
        elif gene in ("uidA", "aggR"):
            true_seqs.append((record.id, record.sequence))
    assert pic_id is not None and pic_seq is not None and len(true_seqs) == 2

    true_assembly = tmp_path / "true.fa"
    with true_assembly.open("w", encoding="utf-8") as handle:
        for seqid, sequence in [(pic_id, pic_seq), *true_seqs]:
            handle.write(f">{seqid} true allele\n{sequence}\n")
    screen_reads_tsv(datadir, true_assembly, "ecoli_dec")
    assert (datadir / "ecoli_dec" / "floors.json").is_file()  # materialized with floors

    # True alleles are all present — the floor does not break true pic calls.
    true_table = screen_reads_tsv(datadir, true_assembly, "ecoli_dec")
    assert tsv_genes(true_table) == {"pic", "uidA", "aggR"}

    homolog = tmp_path / "spate_like.fa"
    with homolog.open("w", encoding="utf-8") as handle:
        handle.write(f">{pic_id} SPATE homolog ~87%\n{_mutate(pic_seq, 0.13, seed=42)}\n")
        for seqid, sequence in true_seqs:
            handle.write(f">{seqid} true allele\n{sequence}\n")

    floored_table = screen_reads_tsv(datadir, homolog, "ecoli_dec")
    assert tsv_genes(floored_table) == {"uidA", "aggR"}  # pic dropped by its floor

    (datadir / "ecoli_dec" / "floors.json").unlink()
    stripped_table = screen_reads_tsv(datadir, homolog, "ecoli_dec")
    assert tsv_genes(stripped_table) == {"pic", "uidA", "aggR"}  # the over-call returns
