"""Summary output formats: abricate-parity TSV/CSV, gapit.summary/1 JSON, Markdown."""

from datetime import datetime
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from gapit import __version__
from gapit.formats.json import ToolDocument, utc_timestamp
from gapit.summary import ABSENT, FIELDSEP, SummaryMatrix


class SummaryParamsDocument(BaseModel, frozen=True):
    """The summary parameters in effect (metric is explicit)."""

    metric: Literal["%COVERAGE", "%IDENTITY"]
    nopath: bool


class SummaryRowDocument(BaseModel, frozen=True):
    """One matrix row: label, distinct-gene count, per-gene value lists."""

    file: str
    num_found: int
    cells: dict[str, list[str]]


class SummaryDocument(BaseModel, frozen=True):
    """gapit.summary/1 — machine-readable summary matrix."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.summary/1"] = Field(default="gapit.summary/1", alias="schema")
    tool: ToolDocument = ToolDocument()
    created_at: str
    params: SummaryParamsDocument
    genes: list[str]
    rows: list[SummaryRowDocument]


def format_summary_tsv(matrix: SummaryMatrix, *, csv: bool) -> str:
    """Render the matrix in abricate --summary shape: '#FILE NUM_FOUND <genes>'
    header, one row per input, ';' cells, '.' absent. Line = sep.join + '\\n'."""
    sep = "," if csv else "\t"
    lines = [sep.join(("#FILE", "NUM_FOUND", *matrix.genes))]
    for row in matrix.rows:
        cells = [
            FIELDSEP.join(row.cells[gene]) if gene in row.cells else ABSENT for gene in matrix.genes
        ]
        lines.append(sep.join((row.file, str(row.num_found), *cells)))
    return "".join(line + "\n" for line in lines)


def render_summary_json(matrix: SummaryMatrix, *, now: datetime) -> str:
    """Serialize the matrix as gapit.summary/1 (indented, schema first)."""
    document = SummaryDocument(
        created_at=utc_timestamp(now),
        params=SummaryParamsDocument(metric=matrix.params.metric, nopath=matrix.params.nopath),
        genes=list(matrix.genes),
        rows=[
            SummaryRowDocument(
                file=row.file,
                num_found=row.num_found,
                cells={gene: list(values) for gene, values in row.cells.items()},
            )
            for row in matrix.rows
        ],
    )
    return document.model_dump_json(indent=2, by_alias=True)


def _md_cell(value: str) -> str:
    """Escape Markdown table metacharacters so cells cannot break the table."""
    return value.replace("\\", "\\\\").replace("|", "\\|")


def render_summary_md(matrix: SummaryMatrix, *, now: datetime) -> str:
    """Render the matrix as deterministic Markdown: YAML frontmatter, then a
    File x gene table with escaped cells ('.' = absent, ';' = multi-value)."""
    columns = ("File", "Num found", *matrix.genes)
    lines: list[str] = [
        "---",
        "schema: gapit.summary/1",
        f"tool: gapit {__version__}",
        f"created_at: {utc_timestamp(now)}",
        f"metric: '{matrix.params.metric}'",
        f"nopath: {str(matrix.params.nopath).lower()}",
        f"files: {len(matrix.rows)}",
        f"genes: {len(matrix.genes)}",
        "---",
        "",
        "# gapit summary matrix",
        "",
        "| " + " | ".join(_md_cell(column) for column in columns) + " |",
        "|" + "---|" * len(columns),
    ]
    for row in matrix.rows:
        cells = [
            _md_cell(FIELDSEP.join(row.cells[gene])) if gene in row.cells else ABSENT
            for gene in matrix.genes
        ]
        lines.append("| " + " | ".join((_md_cell(row.file), str(row.num_found), *cells)) + " |")
    return "\n".join(lines) + "\n"
