"""CLI tests for reads mode (gapit.reads/1) plus committed goldens."""

import json
import re
import shutil
from datetime import UTC, datetime
from pathlib import Path

import pytest
from pydantic import TypeAdapter
from typer.testing import CliRunner

from gapit.cli import app
from gapit.errors import ErrorEnvelope
from gapit.fasta import iter_fasta
from gapit.formats.json import ReadsDocument, render_reads_json
from gapit.formats.md import render_reads_markdown
from gapit.reads import ReadsParams, screen_reads

READS_DB_DIR = Path(__file__).parent / "data" / "reads_db"
READS = Path(__file__).parent / "data" / "reads"
GOLDEN = Path(__file__).parent / "golden"
PINNED_NOW = datetime(2026, 9, 17, 12, 0, 0, tzinfo=UTC)

runner = CliRunner()
reads_adapter = TypeAdapter(ReadsDocument)
PARAMS = ReadsParams(db="tinyreads")


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    target = tmp_path / "datadir"
    shutil.copytree(READS_DB_DIR, target)
    return target


def envelope(stderr: str) -> dict[str, str]:
    return json.loads([line for line in stderr.splitlines() if line.strip()][-1])


# On CI, GITHUB_ACTIONS makes typer force terminal styling on help output;
# the styled runs split option tokens, so strip SGR escapes before asserting.
ANSI_STYLE = re.compile(r"\x1b\[[0-9;]*m")


def test_reads_default_json(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --r1 with no --format, When run, Then gapit.reads/1 JSON on
    stdout and Screening/Detected chatter on stderr."""
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app, ["screen", "--r1", "tetx_full.fq", "--db", "tinyreads", "--datadir", str(datadir)]
    )
    assert result.exit_code == 0
    document = reads_adapter.validate_json(result.stdout)
    assert document.schema_name == "gapit.reads/1"
    assert document.files[0].reads == ["tetx_full.fq"]
    (entry,) = document.files[0].genes
    assert entry.gene == "tetX"
    assert entry.present is True
    assert "Screening reads:" in result.stderr
    assert "Detected 1 present genes" in result.stderr


def test_reads_and_positional_files_are_mutually_exclusive(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "tetx_full.fq",
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
            "junk.fq",
        ],
    )
    assert result.exit_code == 2
    assert envelope(result.stderr)["code"] == "USAGE_ERROR"


def test_reads_mode_rejects_tsv_and_csv(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(READS)
    for fmt in ("tsv", "csv"):
        result = runner.invoke(
            app,
            [
                "screen",
                "--r1",
                "tetx_full.fq",
                "--db",
                "tinyreads",
                "--datadir",
                str(datadir),
                "--format",
                fmt,
            ],
        )
        assert result.exit_code == 2
        assert envelope(result.stderr)["code"] == "USAGE_ERROR"


def test_missing_reads_file_exits_5(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        ["screen", "--r1", "nope.fq", "--db", "tinyreads", "--datadir", str(datadir)],
    )
    assert result.exit_code == 5
    assert envelope(result.stderr)["code"] == "INPUT_NOT_FOUND"


def test_reads_md_output(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "tetx_full.fq",
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
            "--format",
            "md",
        ],
    )
    assert result.exit_code == 0
    assert result.stdout.startswith("---\n")
    assert "| Gene | Breadth% | Depth | Reads | Present |" in result.stdout


def test_r2_without_r1_exits_2(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        ["screen", "--r2", "tetx_R2.fq", "--db", "tinyreads", "--datadir", str(datadir)],
    )
    assert result.exit_code == 2
    assert envelope(result.stderr)["code"] == "USAGE_ERROR"


def test_r2_lane_count_mismatch_exits_2(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "tetx_lane1.fq,tetx_lane2.fq",
            "--r2",
            "tetx_R2.fq",
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
        ],
    )
    assert result.exit_code == 2
    assert envelope(result.stderr)["code"] == "USAGE_ERROR"


def test_empty_comma_element_exits_2(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "tetx_lane1.fq,,tetx_lane2.fq",
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
        ],
    )
    assert result.exit_code == 2
    assert envelope(result.stderr)["code"] == "USAGE_ERROR"


def test_multi_lane_se_union(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given two comma-listed single-end lanes that are each sub-threshold,
    When screened as one sample, Then one gene entry whose union crosses the
    threshold and whose reads_mapped counts both lanes."""
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "tetx_lane1.fq,tetx_lane2.fq",
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
        ],
    )
    assert result.exit_code == 0
    document = reads_adapter.validate_json(result.stdout)
    (entry,) = document.files[0].genes
    assert entry.gene == "tetX"
    assert entry.present is True
    assert entry.reads_mapped == 12
    assert document.files[0].reads == ["tetx_lane1.fq", "tetx_lane2.fq"]


def test_multi_lane_pe_union(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given two paired lanes (2x2 files), When screened, Then one gene entry
    found once with all four files' reads counted."""
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "tetx_pe1_R1.fq,tetx_pe2_R1.fq",
            "--r2",
            "tetx_pe1_R2.fq,tetx_pe2_R2.fq",
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
        ],
    )
    assert result.exit_code == 0
    document = reads_adapter.validate_json(result.stdout)
    (entry,) = document.files[0].genes
    assert entry.present is True
    assert entry.reads_mapped == 12
    assert document.files[0].reads == [
        "tetx_pe1_R1.fq",
        "tetx_pe2_R1.fq",
        "tetx_pe1_R2.fq",
        "tetx_pe2_R2.fq",
    ]


def test_reads_flag_is_gone_from_help() -> None:
    """Given --help for screen, When inspected, Then --reads is absent and
    --r1/--r2 are present."""
    result = runner.invoke(app, ["screen", "--help"], env={"COLUMNS": "100"})
    assert result.exit_code == 0
    help_text = ANSI_STYLE.sub("", result.stdout)
    assert "--reads" not in help_text
    assert "--r1" in help_text
    assert "--r2" in help_text


def test_reads_debug_echoes_minimap2_argv(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --r1 with --debug, When run, Then the minimap2 argv is echoed
    to stderr as a `gapit: run:` line and stdout stays identical to the
    plain run (debug is stderr-only); the plain run emits no argv lines."""
    args = ["screen", "--r1", "tetx_full.fq", "--db", "tinyreads", "--datadir", str(datadir)]
    monkeypatch.chdir(READS)
    debug_run = runner.invoke(app, [*args, "--debug"])
    plain_run = runner.invoke(app, args)
    assert debug_run.exit_code == 0
    run_lines = [line for line in debug_run.stderr.splitlines() if line.startswith("gapit: run:")]
    assert run_lines and run_lines[0].startswith("gapit: run: minimap2 -x sr ")
    assert debug_run.stdout == plain_run.stdout
    assert "gapit: run:" not in plain_run.stderr


def test_golden_reads_json(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Given the tetx fixture run and a pinned now, When rendered, Then the
    JSON is byte-identical to the committed golden."""
    target = tmp_path / "datadir"
    shutil.copytree(READS_DB_DIR, target)
    from gapit.db import discover_databases

    (database,) = discover_databases(target)
    monkeypatch.chdir(READS)
    report = screen_reads(
        [(Path("tetx_full.fq"), None)], database, read_type="sr", min_breadth=90.0, threads=1
    )
    output = render_reads_json([report], PARAMS, now=PINNED_NOW)
    assert output == (GOLDEN / "reads_tinyamr.json").read_text(encoding="utf-8")
    assert reads_adapter.validate_json(output).schema_name == "gapit.reads/1"


def test_golden_reads_md(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    target = tmp_path / "datadir"
    shutil.copytree(READS_DB_DIR, target)
    from gapit.db import discover_databases

    (database,) = discover_databases(target)
    monkeypatch.chdir(READS)
    report = screen_reads(
        [(Path("tetx_full.fq"), None)], database, read_type="sr", min_breadth=90.0, threads=1
    )
    output = render_reads_markdown([report], PARAMS, now=PINNED_NOW)
    assert output == (GOLDEN / "reads_tinyamr.md").read_text(encoding="utf-8")


def _tetx_assembly(tmp_path: Path) -> Path:
    """A single-contig assembly FASTA whose record IS the 522 nt tetX gene
    (the tinyreads sequences file is itself a FASTA of the genes)."""
    tetx = next(
        record
        for record in iter_fasta(READS_DB_DIR / "tinyreads" / "sequences")
        if record.id.startswith("tinyreads~~~tetX")
    )
    assembly = tmp_path / "assembly.fa"
    assembly.write_text(f">contig1\n{tetx.sequence}\n", encoding="utf-8")
    return assembly


def test_fasta_assembly_screens_with_map_ont(datadir: Path, tmp_path: Path) -> None:
    """Given --r1 pointing at an assembly FASTA (the tetX gene as a contig)
    with no --read-type, When screened, Then minimap2 runs with -x map-ont,
    gapit.reads/1 reports params.read_type map-ont with the gene present, and
    the detection note lands on stderr only."""
    assembly = _tetx_assembly(tmp_path)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            str(assembly),
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
            "--debug",
        ],
    )
    assert result.exit_code == 0
    run_lines = [line for line in result.stderr.splitlines() if line.startswith("gapit: run:")]
    assert run_lines and run_lines[0].startswith("gapit: run: minimap2 -x map-ont ")
    document = reads_adapter.validate_json(result.stdout)
    assert document.schema_name == "gapit.reads/1"
    assert document.params.read_type == "map-ont"
    (entry,) = document.files[0].genes
    assert entry.gene == "tetX"
    assert entry.present is True
    assert entry.breadth_pct >= 95.0
    assert entry.reads_mapped == 1
    assert "assembly FASTA detected; using map-ont" in result.stderr
    assert "assembly FASTA" not in result.stdout


def test_fasta_assembly_note_suppressed_by_quiet(datadir: Path, tmp_path: Path) -> None:
    """Given the same assembly run with and without --quiet, When compared,
    Then the detection note is stderr-only and vanishes under --quiet while
    stdout stays byte-identical."""
    assembly = _tetx_assembly(tmp_path)
    args = ["screen", "--r1", str(assembly), "--db", "tinyreads", "--datadir", str(datadir)]
    plain = runner.invoke(app, args)
    quiet = runner.invoke(app, [*args, "--quiet"])
    assert plain.exit_code == 0
    assert quiet.exit_code == 0
    assert "assembly FASTA detected; using map-ont" in plain.stderr
    assert "assembly FASTA detected" not in quiet.stderr
    assert quiet.stdout == plain.stdout


def test_fasta_assembly_explicit_map_ont_is_silent(datadir: Path, tmp_path: Path) -> None:
    """Given --read-type map-ont explicit on an assembly FASTA, When run,
    Then it succeeds with no detection note (the user already chose)."""
    assembly = _tetx_assembly(tmp_path)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            str(assembly),
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
            "--read-type",
            "map-ont",
        ],
    )
    assert result.exit_code == 0
    assert "assembly FASTA detected" not in result.stderr
    document = reads_adapter.validate_json(result.stdout)
    assert document.params.read_type == "map-ont"


def test_fasta_assembly_explicit_sr_exits_2(datadir: Path, tmp_path: Path) -> None:
    """Given an assembly FASTA with explicit --read-type sr, When run, Then
    usage error exit 2 with the frozen message."""
    assembly = _tetx_assembly(tmp_path)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            str(assembly),
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
            "--read-type",
            "sr",
        ],
    )
    assert result.exit_code == 2
    error = envelope(result.stderr)
    assert error["code"] == "USAGE_ERROR"
    assert error["message"] == "assembly FASTA requires map-ont"


def test_fasta_assembly_explicit_map_hifi_exits_2(datadir: Path, tmp_path: Path) -> None:
    """Given an assembly FASTA with explicit --read-type map-hifi, When run,
    Then usage error exit 2 (only map-ont is valid for assemblies)."""
    assembly = _tetx_assembly(tmp_path)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            str(assembly),
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
            "--read-type",
            "map-hifi",
        ],
    )
    assert result.exit_code == 2
    error = envelope(result.stderr)
    assert error["code"] == "USAGE_ERROR"
    assert error["message"] == "assembly FASTA requires map-ont"


def test_mixed_fasta_and_fastq_r1_exits_2(datadir: Path, tmp_path: Path) -> None:
    """Given --r1 listing one FASTA and one FASTQ, When run, Then usage error
    exit 2 before any mapping happens."""
    assembly = _tetx_assembly(tmp_path)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            f"{assembly},{READS / 'tetx_full.fq'}",
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
        ],
    )
    assert result.exit_code == 2
    error = envelope(result.stderr)
    assert error["code"] == "USAGE_ERROR"
    assert error["message"] == "mixed FASTA and FASTQ inputs"


def test_fasta_r1_with_r2_exits_2(datadir: Path, tmp_path: Path) -> None:
    """Given a FASTA --r1 together with --r2, When run, Then usage error exit
    2 (paired-end is FASTQ-only)."""
    assembly = _tetx_assembly(tmp_path)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            str(assembly),
            "--r2",
            str(READS / "tetx_R2.fq"),
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
        ],
    )
    assert result.exit_code == 2
    error = envelope(result.stderr)
    assert error["code"] == "USAGE_ERROR"
    assert error["message"] == "paired-end requires FASTQ"


def test_garbage_reads_file_exits_5(datadir: Path, tmp_path: Path) -> None:
    """Given a --r1 file whose first byte is neither '>' nor '@', When run,
    Then input error exit 5 with the typed code and file context."""
    garbage = tmp_path / "garbage.txt"
    garbage.write_text("Nonsense, not sequencing data\n", encoding="utf-8")
    result = runner.invoke(
        app,
        ["screen", "--r1", str(garbage), "--db", "tinyreads", "--datadir", str(datadir)],
    )
    assert result.exit_code == 5
    error = TypeAdapter(ErrorEnvelope).validate_json(result.stderr.splitlines()[-1])
    assert error.code == "INVALID_READS_FORMAT"
    assert error.context["file"] == str(garbage)


def test_fastq_default_preset_stays_sr(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given a FASTQ --r1 with no --read-type (regression), When run, Then
    params.read_type stays sr."""
    monkeypatch.chdir(READS)
    result = runner.invoke(
        app, ["screen", "--r1", "tetx_full.fq", "--db", "tinyreads", "--datadir", str(datadir)]
    )
    assert result.exit_code == 0
    document = reads_adapter.validate_json(result.stdout)
    assert document.params.read_type == "sr"
