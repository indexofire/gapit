"""Screen result table parsing: contig and reads tables → :class:`ScreenRow`.

The `gapit typing` input boundary. Contig tables (gapit screen / abricate)
carry ``FILE``/``GENE``/``%IDENTITY``/``%COVERAGE``; reads tables
(``gapit screen --r1/--r2``) carry ``#SAMPLE``/``GENE``/``BREADTH%`` and,
only under reads/2's opt-in identity filter, ``IDENTITY%``. Both shapes fold
into the same row model: breadth is the coverage metric, and reads/1's
missing identity arrives as 100.0 because that format cannot gate identity —
the reads engine has already applied its own presence rule before rendering.
"""

import csv
import io
from collections.abc import Iterator
from typing import NamedTuple

from gapit.errors import InputError

_FILE_COLUMNS = ("#FILE", "FILE")
_SAMPLE_COLUMNS = ("#SAMPLE", "SAMPLE")
_READS_MARKERS = (*_SAMPLE_COLUMNS, "BREADTH%")
_READS1_IDENTITY_PCT = 100.0


class ScreenRow(NamedTuple):
    """One data row of a screen result table, reduced to the typing keys."""

    file: str
    gene: str
    database: str
    identity_pct: float
    coverage_pct: float


class _Columns(NamedTuple):
    """Header-name → column-index slots the row parser reads.

    ``identity`` is None only for reads/1 (no IDENTITY% column); ``present``
    is the reads table's presence gate (None on contig tables, whose rows are
    hits by definition)."""

    file: int
    gene: int
    database: int
    identity: int | None
    coverage: int
    present: int | None
    reads: bool

    @property
    def last(self) -> int:
        indexes = [self.file, self.gene, self.database, self.coverage]
        indexes.extend(index for index in (self.identity, self.present) if index is not None)
        return max(indexes)


def _malformed(label: str, line_number: int, detail: str) -> InputError:
    return InputError(
        f"malformed screen result row: {detail} ({label}:{line_number})",
        code="SCREEN_TABLE_MALFORMED",
        context={"file": label, "line": str(line_number)},
    )


def _column_at(indexes: dict[str, int], names: tuple[str, ...]) -> int | None:
    """The first matching header column (``#FILE`` or ``FILE``)."""
    return next((indexes[name] for name in names if name in indexes), None)


def _header(label: str, line_number: int, columns: list[str]) -> _Columns:
    """The first line parsed as the column map (missing slots are malformed)."""
    indexes = dict(zip(columns, range(len(columns)), strict=True))
    if any(marker in indexes for marker in _READS_MARKERS):
        sample_at = _column_at(indexes, _SAMPLE_COLUMNS)
        gene_at = indexes.get("GENE")
        database_at = indexes.get("DATABASE")
        breadth_at = indexes.get("BREADTH%")
        missing = [
            name
            for name, at in (
                ("#SAMPLE", sample_at),
                ("GENE", gene_at),
                ("DATABASE", database_at),
                ("BREADTH%", breadth_at),
            )
            if at is None
        ]
        if sample_at is None or gene_at is None or database_at is None or breadth_at is None:
            raise _malformed(label, line_number, f"header lacks {', '.join(missing)} column(s)")
        return _Columns(
            sample_at,
            gene_at,
            database_at,
            indexes.get("IDENTITY%"),
            breadth_at,
            indexes.get("PRESENT"),
            reads=True,
        )
    file_at = _column_at(indexes, _FILE_COLUMNS)
    gene_at = indexes.get("GENE")
    database_at = indexes.get("DATABASE")
    identity_at = indexes.get("%IDENTITY")
    coverage_at = indexes.get("%COVERAGE")
    missing = [
        name
        for name, at in (
            ("FILE", file_at),
            ("GENE", gene_at),
            ("DATABASE", database_at),
            ("%IDENTITY", identity_at),
            ("%COVERAGE", coverage_at),
        )
        if at is None
    ]
    if (
        file_at is None
        or gene_at is None
        or database_at is None
        or identity_at is None
        or coverage_at is None
    ):
        raise _malformed(label, line_number, f"header lacks {', '.join(missing)} column(s)")
    return _Columns(file_at, gene_at, database_at, identity_at, coverage_at, None, reads=False)


def _iter_rows(text: str, separator: str) -> Iterator[tuple[int, list[str]]]:
    """Physical line number + cells: csv.reader for comma tables, raw split
    for tab tables (the abricate-parity convention)."""
    if separator == ",":
        reader = csv.reader(io.StringIO(text))
        for columns in reader:
            if columns:
                yield reader.line_num, columns
        return
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    for line_number, line in enumerate(lines, start=1):
        yield line_number, line.split(separator)


def _reads_present(cell: str, label: str, line_number: int) -> bool:
    """The reads table's PRESENT cell (yes/no); anything else is malformed."""
    value = cell.strip().lower()
    if value == "yes":
        return True
    if value == "no":
        return False
    raise _malformed(label, line_number, f"PRESENT must be yes or no, got {cell!r}")


def _numeric_detail(header: _Columns) -> str:
    if header.reads:
        return (
            "non-numeric BREADTH%/%IDENTITY"
            if header.identity is not None
            else "non-numeric BREADTH%"
        )
    return "non-numeric %IDENTITY/%COVERAGE"


def parse_table(label: str, text: str) -> list[ScreenRow]:
    """One table's data rows (``label`` names it in errors); the first line
    is the column map and later ``#``-prefixed lines are skipped headers
    (multi-run concatenations). Reads tables skip PRESENT=no rows, so
    ``--all-genes`` output types exactly its present calls."""
    first = text.split("\n", 1)[0]
    separator = "\t" if "\t" in first else "," if "," in first else "\t"
    header: _Columns | None = None
    rows: list[ScreenRow] = []
    for line_number, columns in _iter_rows(text, separator):
        if header is None:
            header = _header(label, line_number, columns)
            continue
        if columns[0].startswith("#"):
            continue
        if len(columns) <= header.last:
            raise _malformed(label, line_number, "expected more columns than present")
        if header.present is not None and not _reads_present(
            columns[header.present], label, line_number
        ):
            continue
        try:
            coverage = float(columns[header.coverage])
            identity = (
                _READS1_IDENTITY_PCT if header.identity is None else float(columns[header.identity])
            )
        except ValueError as exc:
            raise _malformed(label, line_number, _numeric_detail(header)) from exc
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
