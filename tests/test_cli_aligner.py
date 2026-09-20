"""CLI tests for the --aligner engine selector (frozen matrix, 2026-09-19)."""

import json
import re
import shutil
from pathlib import Path

import pytest
from pydantic import TypeAdapter
from typer.testing import CliRunner, Result

from gapit.cli import app
from gapit.db import make_blast_db
from gapit.formats.json import ReadsDocument

AMR_DB_DIR = Path(__file__).parent / "data" / "db"
CONTIGS = Path(__file__).parent / "data" / "contigs"
READS_DB_DIR = Path(__file__).parent / "data" / "reads_db"
READS = Path(__file__).parent / "data" / "reads"
ASSEMBLY = READS_DB_DIR / "tinyreads" / "sequences"

runner = CliRunner()
reads_adapter = TypeAdapter(ReadsDocument)

# On CI, GITHUB_ACTIONS makes typer force terminal styling on help output;
# the styled runs split option tokens, so strip SGR escapes before asserting.
ANSI_STYLE = re.compile(r"\x1b\[[0-9;]*m")


@pytest.fixture()
def amr_datadir(tmp_path: Path) -> Path:
    """Fresh tinyamr datadir with a BLAST index (blastn engine rows)."""
    target = tmp_path / "amr"
    shutil.copytree(AMR_DB_DIR, target)
    make_blast_db(target / "tinyamr" / "sequences", "tinyamr")
    return target


@pytest.fixture()
def reads_datadir(tmp_path: Path) -> Path:
    """Fresh tinyreads datadir (minimap2 engine rows need no BLAST index)."""
    target = tmp_path / "reads"
    shutil.copytree(READS_DB_DIR, target)
    return target


def screen_amr(datadir: Path, *extra: str) -> Result:
    return runner.invoke(app, ["screen", "--db", "tinyamr", "--datadir", str(datadir), *extra])


def screen_reads(datadir: Path, *extra: str) -> Result:
    return runner.invoke(app, ["screen", "--db", "tinyreads", "--datadir", str(datadir), *extra])


def envelope(stderr: str) -> dict[str, str]:
    return json.loads([line for line in stderr.splitlines() if line.strip()][-1])


def reads_document(result: Result) -> dict[str, object]:
    """Parsed gapit.reads/1 document with the wall-clock field dropped."""
    document = reads_adapter.validate_json(result.stdout).model_dump()
    assert isinstance(document, dict)
    document.pop("created_at", None)
    return document


def test_matrix_positional_files_default_to_blastn(amr_datadir: Path) -> None:
    """Given positional contig files with no --aligner (row 1), When screened,
    Then the blastn contig pipeline runs unchanged (abricate TSV on stdout)."""
    result = screen_amr(amr_datadir, "--nopath", str(CONTIGS / "full.fa"))
    assert result.exit_code == 0
    lines = result.stdout.splitlines()
    assert lines[0].startswith("#FILE\t")
    assert lines[1].startswith("full.fa\tcontig1\t1\t79\t+\ttetA\t")


def test_matrix_positional_fasta_with_minimap2_runs_reads_engine(
    reads_datadir: Path,
) -> None:
    """Given a positional assembly FASTA with --aligner minimap2 (row 2), When
    screened, Then the minimap2 engine runs: gapit.reads/1 JSON with the
    resolved map-ont preset and the detection note on stderr only."""
    result = screen_reads(reads_datadir, "--aligner", "minimap2", str(ASSEMBLY))
    assert result.exit_code == 0
    document = reads_adapter.validate_json(result.stdout)
    assert document.schema_name == "gapit.reads/1"
    assert document.params.read_type == "map-ont"
    assert document.files[0].reads == [str(ASSEMBLY)]
    genes = {entry.gene: entry.present for entry in document.files[0].genes}
    assert genes == {"tetX": True, "sulY": True}
    assert "assembly FASTA detected; using map-ont" in result.stderr
    assert "assembly FASTA" not in result.stdout


def test_matrix_positional_fastq_with_minimap2_exits_2(reads_datadir: Path) -> None:
    """Given a positional FASTQ file with --aligner minimap2 (row 2, non-FASTA
    content), When screened, Then usage error exit 2 with the frozen message."""
    result = screen_reads(reads_datadir, "--aligner", "minimap2", str(READS / "tetx_full.fq"))
    assert result.exit_code == 2
    error = envelope(result.stderr)
    assert error["code"] == "USAGE_ERROR"
    assert error["message"] == "minimap2 engine requires FASTA assemblies"


def test_matrix_positional_garbage_with_minimap2_exits_5(
    reads_datadir: Path, tmp_path: Path
) -> None:
    """Given a positional file that is neither FASTA nor FASTQ with --aligner
    minimap2, When screened, Then the typed input error fires (exit 5), not a
    usage error: content detection rejects it before the engine check."""
    garbage = tmp_path / "junk.txt"
    garbage.write_text("not sequencing data\n", encoding="utf-8")
    result = screen_reads(reads_datadir, "--aligner", "minimap2", str(garbage))
    assert result.exit_code == 5
    assert envelope(result.stderr)["code"] == "INVALID_READS_FORMAT"


def test_matrix_positional_blastn_explicit_matches_default(amr_datadir: Path) -> None:
    """Given positional files with explicit --aligner blastn (row 3), When
    compared to the default run, Then stdout is byte-identical (the explicit
    default is a no-op; gbk/embl keep flowing through native normalization)."""
    explicit = screen_amr(amr_datadir, "--aligner", "blastn", "--nopath", str(CONTIGS / "full.fa"))
    default = screen_amr(amr_datadir, "--nopath", str(CONTIGS / "full.fa"))
    assert explicit.exit_code == 0
    assert default.exit_code == 0
    assert explicit.stdout == default.stdout


def test_matrix_reads_default_is_minimap2(
    reads_datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given --r1 FASTQ with no --aligner (row 4), When screened, Then the
    reads path runs unchanged (gapit.reads/1, sr preset)."""
    monkeypatch.chdir(READS)
    result = screen_reads(reads_datadir, "--r1", "tetx_full.fq")
    assert result.exit_code == 0
    document = reads_adapter.validate_json(result.stdout)
    assert document.schema_name == "gapit.reads/1"
    assert document.params.read_type == "sr"


def test_matrix_reads_minimap2_explicit_is_equivalent(
    reads_datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given --r1 with explicit --aligner minimap2 (row 5, redundant but
    valid), When compared to the default run, Then the reads documents are
    identical apart from the wall-clock timestamp."""
    monkeypatch.chdir(READS)
    explicit = screen_reads(reads_datadir, "--r1", "tetx_full.fq", "--aligner", "minimap2")
    default = screen_reads(reads_datadir, "--r1", "tetx_full.fq")
    assert explicit.exit_code == 0
    assert default.exit_code == 0
    assert reads_document(explicit) == reads_document(default)


def test_matrix_reads_blastn_exits_2(reads_datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --r1 with --aligner blastn (row 6), When screened, Then usage
    error exit 2 with the frozen message."""
    monkeypatch.chdir(READS)
    result = screen_reads(reads_datadir, "--r1", "tetx_full.fq", "--aligner", "blastn")
    assert result.exit_code == 2
    error = envelope(result.stderr)
    assert error["code"] == "USAGE_ERROR"
    assert error["message"] == "--aligner blastn is not available for --r1/--r2 reads input"


def test_minimap2_positional_read_type_sr_exits_2(reads_datadir: Path) -> None:
    """Given --aligner minimap2 with explicit --read-type sr on a positional
    FASTA, When screened, Then the frozen assembly-preset rule fires."""
    result = screen_reads(
        reads_datadir, "--aligner", "minimap2", "--read-type", "sr", str(ASSEMBLY)
    )
    assert result.exit_code == 2
    error = envelope(result.stderr)
    assert error["code"] == "USAGE_ERROR"
    assert error["message"] == "assembly FASTA requires map-ont"


def test_minimap2_positional_read_type_map_ont_is_silent(reads_datadir: Path) -> None:
    """Given --aligner minimap2 with explicit --read-type map-ont on a
    positional FASTA, When screened, Then it succeeds with no detection note
    (the user already chose the preset)."""
    result = screen_reads(
        reads_datadir, "--aligner", "minimap2", "--read-type", "map-ont", str(ASSEMBLY)
    )
    assert result.exit_code == 0
    assert "assembly FASTA detected" not in result.stderr
    document = reads_adapter.validate_json(result.stdout)
    assert document.params.read_type == "map-ont"


def test_minimap2_positional_min_breadth_applies(reads_datadir: Path) -> None:
    """Given --aligner minimap2 with --min-breadth 96 on the fixture assembly
    (tetX 97.79% breadth, sulY 95.79%), When screened, Then only tetX is
    present: --min-breadth flows through the positional route."""
    result = screen_reads(
        reads_datadir, "--aligner", "minimap2", "--min-breadth", "96", str(ASSEMBLY)
    )
    assert result.exit_code == 0
    document = reads_adapter.validate_json(result.stdout)
    genes = {entry.gene: entry.present for entry in document.files[0].genes}
    assert genes == {"tetX": True, "sulY": False}


def test_minimap2_positional_format_tsv_exits_2(reads_datadir: Path) -> None:
    """Given --aligner minimap2 with --format tsv, When screened, Then the
    reads output contract rejects it (json is the default)."""
    result = screen_reads(reads_datadir, "--aligner", "minimap2", "--format", "tsv", str(ASSEMBLY))
    assert result.exit_code == 2
    error = envelope(result.stderr)
    assert error["code"] == "USAGE_ERROR"
    assert error["message"] == "--format tsv|csv is not available in reads mode (use json or md)"


def test_minimap2_positional_missing_file_exits_5(reads_datadir: Path) -> None:
    """Given --aligner minimap2 with a missing positional file, When screened,
    Then the typed input error fires before content detection (exit 5)."""
    result = screen_reads(reads_datadir, "--aligner", "minimap2", "nope.fa")
    assert result.exit_code == 5
    assert envelope(result.stderr)["code"] == "INPUT_NOT_FOUND"


def test_minimap2_positional_fofn_exits_2(reads_datadir: Path, tmp_path: Path) -> None:
    """Given --aligner minimap2 together with --fofn, When screened, Then
    usage error exit 2: --fofn is a blastn-engine-only input source and is
    rejected instead of silently ignored."""
    fofn = tmp_path / "files.txt"
    fofn.write_text(f"{ASSEMBLY}\n", encoding="utf-8")
    result = screen_reads(reads_datadir, "--aligner", "minimap2", "--fofn", str(fofn))
    assert result.exit_code == 2
    error = envelope(result.stderr)
    assert error["code"] == "USAGE_ERROR"
    assert error["message"] == "--fofn is not available with --aligner minimap2"


def test_minimap2_without_files_exits_2(reads_datadir: Path) -> None:
    """Given --aligner minimap2 with no positional files, When screened, Then
    usage error exit 2 asking for positional FILEs."""
    result = screen_reads(reads_datadir, "--aligner", "minimap2")
    assert result.exit_code == 2
    error = envelope(result.stderr)
    assert error["code"] == "USAGE_ERROR"
    assert error["message"] == "no input files given (positional FILEs)"


def test_help_lists_aligner_and_no_csv() -> None:
    """Given screen --help, When inspected, Then --aligner is listed with the
    frozen help text and the removed --csv flag is gone."""
    result = runner.invoke(app, ["screen", "--help"], env={"COLUMNS": "100"})
    assert result.exit_code == 0
    help_text = ANSI_STYLE.sub("", result.stdout)
    assert "--aligner" in help_text
    assert "Alignment engine" in help_text
    assert "--csv" not in help_text
