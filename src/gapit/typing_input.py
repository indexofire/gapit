"""`gapit typing` input side: screen report tables → per-FILE gene calls.

Reads gapit/abricate screen TSV (or CSV, auto-detected per file) result
tables and folds them into the evaluation input: one present
:class:`gapit.cluster.GeneCall` per (FILE, GENE) — the best row by
(%IDENTITY, %COVERAGE) descending, the first row winning ties, mirroring
the per-gene best-hit fold the inline engine applied to report hits
(gapit.typing_gene). Rows key on their FILE column, so one designation can
span several tables. Tables arrive as paths or as already-read text (the
piped-stdin table the CLI hands in, labeled ``-``). Validation is typed at
this boundary: a single DATABASE value across every row of every file
(mixed values are DATABASE_MISMATCH), at least one data row
(TYPING_NO_DATA — the database to evaluate is only knowable from rows),
and well-formed numeric columns (SCREEN_TABLE_MALFORMED).
"""

from pathlib import Path
from typing import NamedTuple

from pydantic import BaseModel

from gapit.cluster import GeneCall
from gapit.errors import InputError

_FILE_COLUMNS = ("#FILE", "FILE")

#: The table label for stdin input — the conventional ``-`` marker, also
#: how stdin renders in the JSON ``source`` array and error listings.
STDIN_LABEL = "-"


class ScreenRow(NamedTuple):
    """One data row of a screen result table, reduced to the typing keys."""

    file: str
    gene: str
    database: str
    identity_pct: float
    coverage_pct: float


class FileGeneCalls(BaseModel, frozen=True):
    """One screened FILE and its per-gene folded calls (fold order)."""

    file: str
    calls: dict[str, GeneCall]


class TypingInput(BaseModel, frozen=True):
    """The parsed designation input: the database every row screened
    against plus the per-FILE gene calls, in first-appearance order."""

    database: str
    files: tuple[FileGeneCalls, ...]


def _read_table(path: Path) -> str:
    """One result table's text; missing/unreadable/non-UTF-8 is an InputError."""
    if not path.is_file():
        raise InputError(
            f"screen result file not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        )
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        raise InputError(
            f"screen result is not valid UTF-8: {path}",
            code="SCREEN_TABLE_MALFORMED",
            context={"file": str(path)},
        ) from exc


def _detect_separator(text: str) -> str:
    """Tab when the first line carries one, else comma (summary.py rule)."""
    first = text.split("\n", 1)[0]
    if "\t" in first:
        return "\t"
    if "," in first:
        return ","
    return "\t"


def _malformed(label: str, line_number: int, detail: str) -> InputError:
    return InputError(
        f"malformed screen result row: {detail} ({label}:{line_number})",
        code="SCREEN_TABLE_MALFORMED",
        context={"file": label, "line": str(line_number)},
    )


def _column_at(indexes: dict[str, int], names: tuple[str, ...]) -> int | None:
    """The first matching header column (``#FILE`` or ``FILE``)."""
    return next((indexes[name] for name in names if name in indexes), None)


class _Columns(NamedTuple):
    """Header-name → column-index slots the row parser reads."""

    file: int
    gene: int
    database: int
    identity: int
    coverage: int

    @property
    def last(self) -> int:
        return max(self.file, self.gene, self.database, self.identity, self.coverage)


def _header(label: str, line_number: int, columns: list[str]) -> _Columns:
    """The first line parsed as the column map (missing slots are malformed)."""
    indexes = dict(zip(columns, range(len(columns)), strict=True))
    file_at = _column_at(indexes, _FILE_COLUMNS)
    gene_at = indexes.get("GENE")
    database_at = indexes.get("DATABASE")
    identity_at = indexes.get("%IDENTITY")
    coverage_at = indexes.get("%COVERAGE")
    slots = {
        "FILE": file_at,
        "GENE": gene_at,
        "DATABASE": database_at,
        "%IDENTITY": identity_at,
        "%COVERAGE": coverage_at,
    }
    missing = [name for name, at in slots.items() if at is None]
    if (
        file_at is None
        or gene_at is None
        or database_at is None
        or identity_at is None
        or coverage_at is None
    ):
        raise _malformed(label, line_number, f"header lacks {', '.join(missing)} column(s)")
    return _Columns(file_at, gene_at, database_at, identity_at, coverage_at)


def _parse_rows(label: str, text: str) -> list[ScreenRow]:
    """One table's data rows (``label`` names it in errors); the first line
    is the column map and later ``#``-prefixed lines are skipped headers
    (multi-run concatenations)."""
    separator = _detect_separator(text)
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    header: _Columns | None = None
    rows: list[ScreenRow] = []
    for line_number, line in enumerate(lines, start=1):
        columns = line.split(separator)
        if header is None:
            header = _header(label, line_number, columns)
            continue
        if columns[0].startswith("#"):
            continue
        if len(columns) <= header.last:
            raise _malformed(label, line_number, "expected more columns than present")
        try:
            identity = float(columns[header.identity])
            coverage = float(columns[header.coverage])
        except ValueError as exc:
            raise _malformed(label, line_number, "non-numeric %IDENTITY/%COVERAGE") from exc
        rows.append(
            ScreenRow(
                columns[header.file],
                columns[header.gene],
                columns[header.database],
                identity,
                coverage,
            )
        )
    return rows


def read_typing_input(paths: list[Path], stdin_text: str | None = None) -> TypingInput:
    """Parse and validate the designation input: every table's rows (the
    listed files plus, when ``stdin_text`` is given, the already-read stdin
    table) folded into per-FILE gene calls, with the single DATABASE they
    share."""
    rows: list[ScreenRow] = []
    labels: list[str] = []
    for path in paths:
        labels.append(str(path))
        rows.extend(_parse_rows(str(path), _read_table(path)))
    if stdin_text is not None:
        labels.append(STDIN_LABEL)
        rows.extend(_parse_rows(STDIN_LABEL, stdin_text))
    if not rows:
        listing = ", ".join(labels)
        raise InputError(
            f"no screen result rows to type in: {listing}",
            code="TYPING_NO_DATA",
            context={"files": listing},
        )
    databases: list[str] = []
    for row in rows:
        if row.database not in databases:
            databases.append(row.database)
    if len(databases) > 1:
        raise InputError(
            "screen results mix databases; type each database's results separately:"
            f" {', '.join(databases)}",
            code="DATABASE_MISMATCH",
            context={"databases": ", ".join(databases)},
        )
    best: dict[str, dict[str, ScreenRow]] = {}
    for row in rows:
        current = best.setdefault(row.file, {}).get(row.gene)
        if current is None or (row.identity_pct, row.coverage_pct) > (
            current.identity_pct,
            current.coverage_pct,
        ):
            best[row.file][row.gene] = row
    return TypingInput(
        database=databases[0],
        files=tuple(
            FileGeneCalls(
                file=file,
                calls={
                    gene: GeneCall(
                        gene_id=gene,
                        start=0,
                        end=0,
                        strand="+",
                        coverage_pct=row.coverage_pct,
                        identity_pct=row.identity_pct,
                        verdict="present",
                    )
                    for gene, row in per_gene.items()
                },
            )
            for file, per_gene in best.items()
        ),
    )
