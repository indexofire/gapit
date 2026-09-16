"""Tests for the gapit.report/1 JSON output."""

import shutil
from datetime import UTC, datetime
from pathlib import Path

import pytest
from pydantic import TypeAdapter

from gapit.blast import BlastRow, screen_file
from gapit.db import make_blast_db
from gapit.formats.json import HitDocument, ReportDocument, render_json
from gapit.hits import process_rows
from gapit.report import Report, ScreeningParams

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"
CONTIGS = Path(__file__).parent / "data" / "contigs"
GOLDEN = Path(__file__).parent / "golden"
MULTI_FILES = ("full.fa", "gap.fa", "none.fa", "sort.fa")
PINNED_NOW = datetime(2026, 9, 17, 12, 0, 0, tzinfo=UTC)
PARAMS = ScreeningParams(db="tinyamr")
report_adapter = TypeAdapter(ReportDocument)


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    target = tmp_path / "datadir"
    shutil.copytree(FIXTURE_DB_DIR, target)
    make_blast_db(target / "tinyamr" / "sequences", "tinyamr")
    return target


@pytest.fixture()
def reports(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> list[Report]:
    """Real-pipeline reports for the four fixture files, run from CONTIGS so
    Report.file is a stable bare filename."""
    monkeypatch.chdir(CONTIGS)
    from gapit.db import Database

    database = Database(
        name="tinyamr", path=datadir / "tinyamr", sequences_path=datadir / "tinyamr" / "sequences"
    )
    return [screen_file(Path(name), database, PARAMS) for name in MULTI_FILES]


def test_golden_json(reports: list[Report]) -> None:
    """Given the multi-file fixture run and a pinned now, When rendered, Then
    the output is byte-identical to the committed golden JSON."""
    assert render_json(reports, PARAMS, now=PINNED_NOW) == (GOLDEN / "screen_multi.json").read_text(
        encoding="utf-8"
    )


def test_schema_literal_comes_first(reports: list[Report]) -> None:
    output = render_json(reports, PARAMS, now=PINNED_NOW)
    assert output.startswith('{\n  "schema": "gapit.report/1",')
    document = report_adapter.validate_json(output)
    assert document.schema_name == "gapit.report/1"


def test_round_trip_parses_into_report_document(reports: list[Report]) -> None:
    """Given emitted JSON, When parsed into ReportDocument, Then it validates."""
    document = report_adapter.validate_json(render_json(reports, PARAMS, now=PINNED_NOW))
    assert document.tool.name == "gapit"
    assert document.created_at == "2026-09-17T12:00:00Z"
    assert document.params.db == "tinyamr"
    assert [f.file for f in document.files] == list(MULTI_FILES)


def test_hit_key_order_is_documented(reports: list[Report]) -> None:
    output = render_json(reports, PARAMS, now=PINNED_NOW)
    keys = list(HitDocument.model_fields)
    assert keys == [
        "sequence",
        "start",
        "end",
        "strand",
        "gene",
        "coverage",
        "coverage_map",
        "gaps",
        "coverage_pct",
        "identity_pct",
        "database",
        "accession",
        "product",
        "resistance",
    ]
    positions = [output.index(f'"{key}"') for key in keys]
    assert positions == sorted(positions)


def test_zero_hit_file_present_with_empty_list(reports: list[Report]) -> None:
    document = report_adapter.validate_json(render_json(reports, PARAMS, now=PINNED_NOW))
    by_file = {f.file: f.hits for f in document.files}
    assert by_file["none.fa"] == []
    assert len(by_file["sort.fa"]) == 4


def test_pct_fields_are_rounded_numbers() -> None:
    """Given coverage_pct 79.996 and identity 96.907, When rendered, Then the
    JSON numbers are rounded to 2 decimals (80.0 / 96.91)."""
    row = BlastRow(
        qseqid="c",
        qstart=1,
        qend=20000,
        qlen=20000,
        sseqid="db~~~g~~~a~~~r",
        sstart=1,
        send=25000,
        slen=25000,
        sstrand="plus",
        evalue=1e-40,
        length=20000,
        pident=96.907,
        gaps=1,
        gapopen=1,
        stitle="db~~~g~~~a~~~r product",
    )
    (hit,) = process_rows([row], mincov=0.0, default_db="db")
    document = report_adapter.validate_json(
        render_json([Report(file="x.fa", hits=(hit,))], ScreeningParams(db="db"), now=PINNED_NOW)
    )
    entry = document.files[0].hits[0]
    assert entry.coverage_pct == 80.0
    assert entry.identity_pct == 96.91


def test_cli_format_json(
    reports: list[Report], datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given --format json on the CLI, When run, Then stdout parses as
    gapit.report/1 (created_at unpinned, just structurally checked)."""
    from typer.testing import CliRunner

    from gapit.cli import app

    monkeypatch.chdir(CONTIGS)
    result = CliRunner().invoke(
        app,
        ["screen", *MULTI_FILES, "--db", "tinyamr", "--datadir", str(datadir), "--format", "json"],
    )
    assert result.exit_code == 0
    document = report_adapter.validate_json(result.stdout)
    assert document.schema_name == "gapit.report/1"
    assert document.created_at.endswith("Z")
