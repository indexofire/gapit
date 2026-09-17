"""Summary mode core: parse abricate-format report tables into a gene matrix.

Mirrors abricate 1.4.0 ``summary_table`` (SPEC.md §6) exactly where it is
defined, and replaces its silent-undef edges with typed InputErrors
(``SUMMARY_MALFORMED``) — a documented [gapit-extension] divergence.
"""

from collections.abc import Callable
from pathlib import Path, PurePath
from typing import Literal

from pydantic import BaseModel

from gapit.errors import InputError

FIELDSEP = ";"
ABSENT = "."

Warn = Callable[[str], None]


class SummaryParams(BaseModel, frozen=True):
    """Summary parameters in effect (metric + path display)."""

    identity: bool = False
    nopath: bool = False

    @property
    def metric(self) -> Literal["%COVERAGE", "%IDENTITY"]:
        """The report column summarized into cells."""
        return "%IDENTITY" if self.identity else "%COVERAGE"


class SummaryRow(BaseModel, frozen=True):
    """One matrix row: a display label, its distinct-gene count, and the
    per-gene cell values as the original report strings in row order."""

    file: str
    num_found: int
    cells: dict[str, tuple[str, ...]]


class SummaryMatrix(BaseModel, frozen=True):
    """Canonical in-memory summary result: sorted gene universe + ordered rows."""

    params: SummaryParams
    genes: tuple[str, ...]
    rows: tuple[SummaryRow, ...]


def _read_report(path: Path) -> str:
    """Read one report table; missing/unreadable/non-UTF-8 is an InputError."""
    if not path.is_file():
        raise InputError(
            f"report file not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        )
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        raise InputError(
            f"report is not valid UTF-8: {path}",
            code="SUMMARY_MALFORMED",
            context={"file": str(path)},
        ) from exc


def _lines(text: str) -> list[str]:
    """Split like Perl's ``while (<$fh>)``: on \\n, no phantom final line."""
    if not text:
        return []
    parts = text.split("\n")
    if parts[-1] == "":
        parts.pop()
    return parts


def _detect_separator(text: str) -> str:
    """[gapit-extension] auto-detect the table separator per file: tab if the
    first line has any, else comma, else tab (single-column degenerate)."""
    first = _lines(text)[0] if text else ""
    if "\t" in first:
        return "\t"
    if "," in first:
        return ","
    return "\t"


def _malformed(path: Path, line_number: int, detail: str) -> InputError:
    return InputError(
        f"malformed report row: {detail} ({path}:{line_number})",
        code="SUMMARY_MALFORMED",
        context={"file": str(path), "line": str(line_number)},
    )


def build_summary(paths: list[Path], params: SummaryParams, *, warn: Warn) -> SummaryMatrix:
    """Aggregate report table(s) into the summary matrix.

    Dutch mode (exactly one input) keys rows by each table row's FILE column;
    otherwise rows are keyed by input filename as given. Rows are sorted by
    that key — not by display label — matching upstream. The first row seen
    anywhere is the column map (upstream ``@hdr``); later '#' rows are skipped.
    """
    dutch = len(paths) == 1
    data: dict[str, dict[str, list[str]]] = {}
    seen: set[str] = set()
    indexes: dict[str, int] | None = None
    for path in paths:
        key = str(path)
        if key in seen:
            warn(f"Skipping duplicate file: {key}")
            continue
        seen.add(key)
        text = _read_report(path)
        if not dutch:
            data[key] = {}
        separator = _detect_separator(text)
        for line_number, line in enumerate(_lines(text), start=1):
            columns = line.split(separator)
            if indexes is None:
                # Header-name -> column index with Perl zip semantics (later dups win).
                indexes = dict(zip(columns, range(len(columns)), strict=True))
            if columns[0].startswith("#"):
                continue
            assert indexes is not None  # set on the first line, before any data row
            gene_at = indexes.get("GENE")
            metric_at = indexes.get(params.metric)
            if gene_at is None:
                raise _malformed(path, line_number, "header has no GENE column")
            if metric_at is None:
                raise _malformed(path, line_number, f"header has no {params.metric} column")
            if len(columns) <= max(gene_at, metric_at):
                raise _malformed(
                    path, line_number, f"expected >= {max(gene_at, metric_at) + 1} columns"
                )
            gene, value = columns[gene_at], columns[metric_at]
            file_key = PurePath(columns[0]).name if params.nopath else columns[0]
            row_key = file_key if dutch else key
            data.setdefault(row_key, {}).setdefault(gene, []).append(value)
    genes = tuple(sorted({gene for hits in data.values() for gene in hits}))
    rows = tuple(
        SummaryRow(
            file=PurePath(row_key).name if params.nopath else row_key,
            num_found=len(data[row_key]),
            cells={gene: tuple(data[row_key][gene]) for gene in genes if gene in data[row_key]},
        )
        for row_key in sorted(data)
    )
    return SummaryMatrix(params=params, genes=genes, rows=rows)
