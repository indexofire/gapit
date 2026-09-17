"""Unit tests for summary output formats: TSV/CSV parity + gapit.summary/1 JSON/MD."""

from datetime import UTC, datetime
from pathlib import Path

import pytest
from pydantic import TypeAdapter

from gapit.formats.summary import (
    SummaryDocument,
    format_summary_tsv,
    render_summary_json,
    render_summary_md,
)
from gapit.summary import SummaryParams, build_summary

FIXTURES = Path(__file__).parent / "data" / "summary"
PINNED_NOW = datetime(2026, 9, 17, 12, 0, 0, tzinfo=UTC)
PARAMS = SummaryParams()
summary_adapter = TypeAdapter(SummaryDocument)

MULTI = ["sample_a.tsv", "sample_b.tsv", "empty.tsv"]

EXPECTED_MULTI_TSV = (
    "#FILE\tNUM_FOUND\tfeature_a\tfeature_b\n"
    "empty.tsv\t0\t.\t.\n"
    "sample_a.tsv\t2\t99.50;52.00\t76.00\n"
    "sample_b.tsv\t2\t90.00\t100.00\n"
)


def build(names: list[str], params: SummaryParams = PARAMS):
    return build_summary([Path(name) for name in names], params, warn=lambda _message: None)


def test_tsv_byte_parity_with_abricate(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given the verified multi-file run, When rendered as TSV, Then output is
    byte-identical to real abricate --summary (incl. ';' cells and '.' absent)."""
    monkeypatch.chdir(FIXTURES)
    assert format_summary_tsv(build(MULTI), csv=False) == EXPECTED_MULTI_TSV


def test_csv_joins_with_commas(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given csv=True, When rendered, Then the same matrix comma-joined
    (upstream --csv output shape)."""
    monkeypatch.chdir(FIXTURES)
    assert format_summary_tsv(build(MULTI), csv=True) == EXPECTED_MULTI_TSV.replace("\t", ",")


def test_tsv_zero_genes_header_has_no_trailing_separator(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Given only empty reports, When rendered, Then the header is
    '#FILE\\tNUM_FOUND' with no trailing tab and no data rows."""
    monkeypatch.chdir(FIXTURES)
    output = format_summary_tsv(build(["empty.tsv"]), csv=False)
    assert output == "#FILE\tNUM_FOUND\n"


def test_json_is_gapit_summary_1(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given a matrix, When rendered as JSON, Then schema/tool/created_at lead,
    params state the metric explicitly, and rows are typed with num_found and
    per-gene value lists preserving the original strings."""
    monkeypatch.chdir(FIXTURES)
    document = summary_adapter.validate_json(render_summary_json(build(MULTI), now=PINNED_NOW))
    assert document.schema_name == "gapit.summary/1"
    assert document.tool.name == "gapit"
    assert document.created_at == "2026-09-17T12:00:00Z"
    assert document.params.metric == "%COVERAGE"
    assert document.params.nopath is False
    assert document.genes == ["feature_a", "feature_b"]
    sample_a = next(row for row in document.rows if row.file == "sample_a.tsv")
    assert sample_a.num_found == 2
    assert sample_a.cells == {"feature_a": ["99.50", "52.00"], "feature_b": ["76.00"]}
    empty = document.rows[0]
    assert empty.file == "empty.tsv"
    assert empty.num_found == 0
    assert empty.cells == {}


def test_json_identity_param_carries_metric(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --identity, When rendered as JSON, Then params.metric is
    %IDENTITY and cells are the identity strings."""
    monkeypatch.chdir(FIXTURES)
    document = summary_adapter.validate_json(
        render_summary_json(build(["sample_a.tsv"], SummaryParams(identity=True)), now=PINNED_NOW)
    )
    assert document.params.metric == "%IDENTITY"
    (row,) = document.rows
    assert row.cells["feature_a"] == ["98.75", "91.00"]


def test_md_frontmatter_and_matrix(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given a matrix, When rendered as Markdown, Then YAML frontmatter carries
    schema/tool/created_at/metric and the table shows one row per file."""
    monkeypatch.chdir(FIXTURES)
    output = render_summary_md(build(MULTI), now=PINNED_NOW)
    assert output.startswith("---\n")
    assert "schema: gapit.summary/1" in output
    assert "tool: gapit" in output
    assert "created_at: 2026-09-17T12:00:00Z" in output
    assert "metric: '%COVERAGE'" in output
    assert "| File | Num found | feature_a | feature_b |" in output
    assert "| sample_a.tsv | 2 | 99.50;52.00 | 76.00 |" in output
    assert "| empty.tsv | 0 | . | . |" in output


def test_md_escapes_pipes_in_cells(tmp_path: Path) -> None:
    """Given a gene label containing '|', When rendered as Markdown, Then the
    pipe is escaped so the table structure survives."""
    header = (FIXTURES / "empty.tsv").read_text(encoding="utf-8")
    text = (
        header + "k.fa\tc\t1\t9\t+\twe|ird\t1-9/9\t===============\t0/0\t90.00\t90.00\tdb\ta\tp\t\n"
    )
    path = tmp_path / "pipe.tsv"
    path.write_text(text, encoding="utf-8")
    output = render_summary_md(build([str(path)]), now=PINNED_NOW)
    assert "we\\|ird" in output
    assert "| we|ird" not in output
