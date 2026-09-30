"""CLI end-to-end tests for --merge-fragments over the fragments fixtures
(real normalization -> blastn pipeline; goldens committed under tests/golden)."""

import json
import shutil
from datetime import UTC, datetime
from pathlib import Path

from typer.testing import CliRunner, Result

from gapit.blast import screen_file
from gapit.cli import app
from gapit.db import Database, make_blast_db
from gapit.formats.json import render_json
from gapit.formats.md import render_markdown
from gapit.report import Report, ScreeningParams

FIXTURES = Path(__file__).parent / "data" / "fragments"
GOLDEN = Path(__file__).parent / "golden"
FILES = ("split4060.fa", "split5050.fa", "split3040.fa")
PINNED_NOW = datetime(2026, 9, 17, 12, 0, 0, tzinfo=UTC)
PARAMS = ScreeningParams(db="fragdb")

runner = CliRunner()


def fragdb(tmp_path: Path) -> Database:
    """Fresh datadir with a built fragdb index (one per test)."""
    datadir = tmp_path / "datadir"
    shutil.copytree(FIXTURES / "fragdb", datadir / "fragdb")
    make_blast_db(datadir / "fragdb" / "sequences", "fragdb")
    return Database(
        name="fragdb", path=datadir / "fragdb", sequences_path=datadir / "fragdb" / "sequences"
    )


def merged_reports(database: Database) -> list[Report]:
    """Merge-mode reports for the three fixtures, bare filenames like the
    goldens (Report.file is the only field rewritten)."""
    reports = [
        screen_file(FIXTURES / name, database, PARAMS, dbtype="nucl", merge_fragments=True)
        for name in FILES
    ]
    return [
        report.model_copy(update={"file": name})
        for report, name in zip(reports, FILES, strict=True)
    ]


def screen(database: Database, *extra: str) -> Result:
    return runner.invoke(
        app,
        [
            "screen",
            "--db",
            "fragdb",
            "--datadir",
            str(database.path.parent),
            *extra,
        ],
    )


def test_golden_tsv_merge_mode(tmp_path: Path) -> None:
    """Given the three fragment fixtures with --merge-fragments --nopath,
    When screened, Then stdout is byte-identical to the committed golden TSV:
    one row per merged gene, none for the 30/40 fail case."""
    database = fragdb(tmp_path)
    result = screen(
        database, "--merge-fragments", "--nopath", *[str(FIXTURES / name) for name in FILES]
    )
    assert result.exit_code == 0
    assert result.stdout == (GOLDEN / "fragments_merge.tsv").read_text(encoding="utf-8")


def test_golden_json_merge_mode(tmp_path: Path) -> None:
    """Given the same run rendered as JSON with a pinned now, When compared,
    Then byte-identical to the golden: merged rows carry ``merged: true`` and
    the additive ``fragments`` array, plain rows carry neither."""
    output = render_json(merged_reports(fragdb(tmp_path)), PARAMS, now=PINNED_NOW)
    assert output == (GOLDEN / "fragments_merge.json").read_text(encoding="utf-8")
    document = json.loads(output)
    hit = document["files"][0]["hits"][0]
    assert hit["merged"] is True
    assert [(f["contig"], f["strand"]) for f in hit["fragments"]] == [
        ("contigA", "+"),
        ("contigB", "-"),
    ]
    assert "merged" not in document["files"][2]["hits"]


def test_golden_md_merge_mode(tmp_path: Path) -> None:
    """Given the same run rendered as Markdown with a pinned now, When
    compared, Then byte-identical to the golden: the table rows plus one
    fragments detail line per merged gene."""
    output = render_markdown(merged_reports(fragdb(tmp_path)), PARAMS, now=PINNED_NOW)
    assert output == (GOLDEN / "fragments_merge.md").read_text(encoding="utf-8")
    assert "- `demoGene` merged from 2 fragments:" in output


def test_default_path_reports_nothing_without_the_flag(tmp_path: Path) -> None:
    """Given the same fixtures WITHOUT --merge-fragments, When screened, Then
    every fragment stays below mincov and the gene is absent — default output
    is byte-identical with and without the opt-in flag spelling."""
    database = fragdb(tmp_path)
    files = [str(FIXTURES / name) for name in FILES]
    default = screen(database, "--nopath", *files)
    off = screen(database, "--no-merge-fragments", "--nopath", *files)
    assert default.exit_code == 0
    assert default.stdout.splitlines() == [default.stdout.splitlines()[0]]
    assert default.stdout == off.stdout


def test_union_below_mincov_stays_unreported_in_merge_mode(tmp_path: Path) -> None:
    """Given the 30/40 fail case (union 70%), When screened with
    --merge-fragments, Then no row: merging only rescues unions that reach
    the threshold."""
    database = fragdb(tmp_path)
    result = screen(database, "--merge-fragments", "--nopath", str(FIXTURES / "split3040.fa"))
    assert result.exit_code == 0
    assert result.stdout.splitlines() == [result.stdout.splitlines()[0]]
    assert "Found 0 genes" in result.stderr


def test_merge_mode_respects_mincov_threshold(tmp_path: Path) -> None:
    """Given the 30/40 case (union exactly 70.0%), When screened with
    --merge-fragments --mincov 70, Then the merged hit appears (>= keeps);
    at --mincov 70.01 it does not."""
    database = fragdb(tmp_path)
    kept = screen(
        database, "--merge-fragments", "--mincov", "70", "--nopath", str(FIXTURES / "split3040.fa")
    )
    dropped = screen(
        database,
        "--merge-fragments",
        "--mincov",
        "70.01",
        "--nopath",
        str(FIXTURES / "split3040.fa"),
    )
    assert len(kept.stdout.splitlines()) == 2
    assert "\t70.00\t" in kept.stdout.splitlines()[1]
    assert len(dropped.stdout.splitlines()) == 1


def test_cli_json_carries_fragments(tmp_path: Path) -> None:
    """Given --merge-fragments --format json on the CLI, When parsed, Then
    the merged hit documents its fragments (contig, spans, strand, pcts)."""
    database = fragdb(tmp_path)
    result = screen(
        database, "--merge-fragments", "--format", "json", str(FIXTURES / "split4060.fa")
    )
    assert result.exit_code == 0
    hit = json.loads(result.stdout)["files"][0]["hits"][0]
    assert hit["merged"] is True
    assert [frag["contig"] for frag in hit["fragments"]] == ["contigA", "contigB"]
    assert hit["fragments"][0]["identity_pct"] == 98.75
    assert hit["identity_pct"] == 99.5


def test_merge_fragments_rejected_in_reads_mode(tmp_path: Path) -> None:
    """Given --merge-fragments with --r1, When screened, Then exit 2 with a
    usage envelope (mirrors the --jobs reads-mode guard)."""
    database = fragdb(tmp_path)
    reads = FIXTURES.parent / "reads"
    result = screen(
        database,
        "--merge-fragments",
        "--r1",
        str(reads / "tetx_R1.fq"),
        "--read-type",
        "sr",
    )
    assert result.exit_code == 2
    assert json.loads(result.stderr.splitlines()[-1])["code"] == "USAGE_ERROR"


def test_merge_fragments_rejected_with_aligner_minimap2(tmp_path: Path) -> None:
    """Given --merge-fragments with --aligner minimap2, When screened, Then
    exit 2 with a usage envelope."""
    database = fragdb(tmp_path)
    result = screen(
        database,
        "--merge-fragments",
        "--aligner",
        "minimap2",
        str(FIXTURES / "split4060.fa"),
    )
    assert result.exit_code == 2
    assert json.loads(result.stderr.splitlines()[-1])["code"] == "USAGE_ERROR"
