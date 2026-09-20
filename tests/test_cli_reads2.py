"""CLI tests for gapit.reads/2: opt-in --min-identity/--min-mapq, document
selection, validation, blastn-engine rejection, schema introspection, and the
reads/1 byte-identity lock when both thresholds are off.
"""

import json
import re
import shutil
from pathlib import Path

import pytest
from pydantic import TypeAdapter
from typer.testing import CliRunner

from gapit.cli import app
from gapit.formats.json import Reads2Document

READS2_DB = Path(__file__).parent / "data" / "reads2_db"
READS2 = Path(__file__).parent / "data" / "reads2"
READS = Path(__file__).parent / "data" / "reads"

runner = CliRunner()
reads2_adapter = TypeAdapter(Reads2Document)

# On CI, GITHUB_ACTIONS makes typer force terminal styling on help output;
# the styled runs split option tokens, so strip SGR escapes before asserting.
ANSI_STYLE = re.compile(r"\x1b\[[0-9;]*m")


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    target = tmp_path / "datadir"
    shutil.copytree(READS2_DB, target)
    return target


def envelope(stderr: str) -> dict[str, str]:
    return json.loads([line for line in stderr.splitlines() if line.strip()][-1])


# --- reads/1 freeze -----------------------------------------------------------


def test_no_flags_emits_reads1(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --r1 with no new flags, When run, Then the document is
    gapit.reads/1 (the frozen default)."""
    monkeypatch.chdir(READS2)
    result = runner.invoke(
        app, ["screen", "--r1", "sr_homologs.fq", "--db", "homologs", "--datadir", str(datadir)]
    )
    assert result.exit_code == 0
    assert json.loads(result.stdout)["schema"] == "gapit.reads/1"


def test_explicit_off_flags_are_byte_identical_to_no_flags(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the same run with no flags vs explicit --min-identity 0
    --min-mapq 0, When compared, Then stdout is byte-identical (off means
    off; the default path gains nothing)."""
    monkeypatch.chdir(READS2)
    args = ["screen", "--r1", "sr_homologs.fq", "--db", "homologs", "--datadir", str(datadir)]
    plain = runner.invoke(app, args)
    explicit_off = runner.invoke(app, [*args, "--min-identity", "0", "--min-mapq", "0"])
    assert plain.exit_code == 0
    assert explicit_off.exit_code == 0
    assert explicit_off.stdout == plain.stdout


def test_off_run_does_not_add_cs_to_minimap2_argv(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the default /1 run, When --debug echoes the minimap2 argv, Then
    --cs is absent (the /1 invocation is unchanged; only /2 mode needs the
    NM tag)."""
    monkeypatch.chdir(READS2)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "sr_homologs.fq",
            "--db",
            "homologs",
            "--datadir",
            str(datadir),
            "--debug",
        ],
    )
    assert result.exit_code == 0
    run_lines = [line for line in result.stderr.splitlines() if line.startswith("gapit: run:")]
    assert run_lines
    assert "--cs" not in run_lines[0]


# --- reads/2 opt-in ------------------------------------------------------------


def test_min_identity_emits_reads2(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --min-identity 95, When run, Then the document is gapit.reads/2
    with the threshold in params and mean_identity_pct on gene entries."""
    monkeypatch.chdir(READS2)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "sr_homologs.fq",
            "--db",
            "homologs",
            "--datadir",
            str(datadir),
            "--min-identity",
            "95",
        ],
    )
    assert result.exit_code == 0
    document = reads2_adapter.validate_json(result.stdout)
    assert document.schema_name == "gapit.reads/2"
    assert document.params.min_identity == 95.0
    assert document.params.min_mapq == 0
    (entry,) = document.files[0].genes
    assert entry.gene == "geneA"
    assert entry.mean_identity_pct == 100.0


def test_min_mapq_emits_reads2(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --min-mapq 30 alone, When run, Then the document is gapit.reads/2
    (either flag opts in) and params carry both thresholds."""
    monkeypatch.chdir(READS2)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "sr_homologs.fq",
            "--db",
            "homologs",
            "--datadir",
            str(datadir),
            "--min-mapq",
            "30",
        ],
    )
    assert result.exit_code == 0
    document = reads2_adapter.validate_json(result.stdout)
    assert document.schema_name == "gapit.reads/2"
    assert document.params.min_mapq == 30
    assert document.params.min_identity == 0.0


def test_min_identity_md_output_is_reads2(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --min-identity with --format md, When run, Then frontmatter says
    gapit.reads/2 and lists the threshold."""
    monkeypatch.chdir(READS2)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "sr_homologs.fq",
            "--db",
            "homologs",
            "--datadir",
            str(datadir),
            "--min-identity",
            "95",
            "--format",
            "md",
        ],
    )
    assert result.exit_code == 0
    assert "schema: gapit.reads/2" in result.stdout
    assert "min_identity: 95.0" in result.stdout
    assert "| Gene | Breadth% |" in result.stdout


def test_min_identity_assembly_route_emits_reads2(datadir: Path, tmp_path: Path) -> None:
    """Given --aligner minimap2 with --min-identity (assembly FASTA), When
    run, Then the reads pipeline applies the filter and emits reads/2."""
    assembly = tmp_path / "assembly.fa"
    first = (READS2_DB / "homologs" / "sequences").read_text(encoding="utf-8").splitlines()
    assembly.write_text(f"{first[0]}\n{first[1]}\n", encoding="utf-8")
    result = runner.invoke(
        app,
        [
            "screen",
            "--aligner",
            "minimap2",
            str(assembly),
            "--db",
            "homologs",
            "--datadir",
            str(datadir),
            "--min-identity",
            "95",
        ],
    )
    assert result.exit_code == 0
    document = reads2_adapter.validate_json(result.stdout)
    assert document.schema_name == "gapit.reads/2"
    assert document.params.read_type == "map-ont"


def test_reads2_debug_argv_contains_cs(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given a /2 run with --debug, When the minimap2 argv is echoed, Then it
    carries --cs (the PAF flag that emits NM:i:)."""
    monkeypatch.chdir(READS2)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "sr_homologs.fq",
            "--db",
            "homologs",
            "--datadir",
            str(datadir),
            "--min-identity",
            "95",
            "--debug",
        ],
    )
    assert result.exit_code == 0
    run_lines = [line for line in result.stderr.splitlines() if line.startswith("gapit: run:")]
    assert run_lines
    assert " --cs " in run_lines[0]


# --- validation + engine rejection ---------------------------------------------


def test_min_identity_above_100_exits_2(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(READS2)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "sr_homologs.fq",
            "--db",
            "homologs",
            "--datadir",
            str(datadir),
            "--min-identity",
            "100.5",
        ],
    )
    assert result.exit_code == 2
    error = envelope(result.stderr)
    assert error["code"] == "USAGE_ERROR"
    assert error["message"] == "--min-identity must be in [0, 100]: got 100.5"


def test_min_identity_negative_exits_2(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(READS2)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "sr_homologs.fq",
            "--db",
            "homologs",
            "--datadir",
            str(datadir),
            "--min-identity",
            "-1",
        ],
    )
    assert result.exit_code == 2
    assert envelope(result.stderr)["code"] == "USAGE_ERROR"


def test_min_mapq_negative_exits_2(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(READS2)
    result = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "sr_homologs.fq",
            "--db",
            "homologs",
            "--datadir",
            str(datadir),
            "--min-mapq",
            "-1",
        ],
    )
    assert result.exit_code == 2
    error = envelope(result.stderr)
    assert error["code"] == "USAGE_ERROR"
    assert error["message"] == "--min-mapq must be >= 0: got -1"


def test_blastn_engine_rejects_identity_flag() -> None:
    """Given positional contig files (blastn engine) with --min-identity,
    When run, Then usage error exit 2 with the frozen reads-mode-only
    message."""
    contigs = Path(__file__).parent / "data" / "contigs" / "full.fa"
    result = runner.invoke(
        app,
        [
            "screen",
            str(contigs),
            "--db",
            "homologs",
            "--datadir",
            str(READS2_DB),
            "--min-identity",
            "90",
        ],
    )
    assert result.exit_code == 2
    error = envelope(result.stderr)
    assert error["code"] == "USAGE_ERROR"
    assert error["message"] == "--min-identity/--min-mapq are reads-mode only (minimap2 engine)"


def test_blastn_engine_rejects_mapq_flag() -> None:
    contigs = Path(__file__).parent / "data" / "contigs" / "full.fa"
    result = runner.invoke(
        app,
        [
            "screen",
            str(contigs),
            "--db",
            "homologs",
            "--datadir",
            str(READS2_DB),
            "--min-mapq",
            "20",
        ],
    )
    assert result.exit_code == 2
    assert envelope(result.stderr)["code"] == "USAGE_ERROR"


# --- schema introspection -------------------------------------------------------


def test_schema_reads2_emits_valid_schema() -> None:
    """Given `gapit schema reads2`, When run, Then valid JSON Schema for the
    gapit.reads/2 document (mean_identity_pct present, /1 const absent)."""
    result = runner.invoke(app, ["schema", "reads2"])
    assert result.exit_code == 0
    schema = json.loads(result.stdout)
    assert schema["title"] == "Reads2Document"
    assert "mean_identity_pct" in result.stdout
    assert "min_identity" in result.stdout


def test_schema_unknown_name_lists_reads2() -> None:
    """Given an unknown schema name, When run, Then the help list of choices
    includes reads2."""
    result = runner.invoke(app, ["schema", "nope"])
    assert result.exit_code == 2
    assert "reads2" in envelope(result.stderr)["message"]


def test_screen_help_lists_new_flags() -> None:
    result = runner.invoke(app, ["screen", "--help"], env={"COLUMNS": "100"})
    assert result.exit_code == 0
    help_text = ANSI_STYLE.sub("", result.stdout)
    assert "--min-identity" in help_text
    assert "--min-mapq" in help_text
