"""`gapit typing` input side: screen report tables → per-FILE gene calls.

Reads gapit/abricate contig screen tables and gapit reads tables (the
``#SAMPLE``/``BREADTH%`` shape), TSV or CSV, auto-detected per file, and
folds them into the evaluation input: one present
:class:`gapit.cluster.GeneCall` per (FILE, GENE) — the best row by
(%IDENTITY, %COVERAGE) descending, the first row winning ties, mirroring
the per-gene best-hit fold the inline engine applied to report hits
(gapit.typing_gene). Rows key on their FILE/#SAMPLE column, so one
designation can span several tables. Tables arrive as paths or as
already-read text (the piped-stdin table the CLI hands in, labeled ``-``).
Validation is typed at this boundary: a single DATABASE value across every
row of every file (mixed values are DATABASE_MISMATCH), at least one data
row (TYPING_NO_DATA — the database to evaluate is only knowable from rows),
and well-formed numeric columns (SCREEN_TABLE_MALFORMED). Table-shape
parsing lives in :mod:`gapit.typing_table`.
"""

from pathlib import Path

from pydantic import BaseModel

from gapit.cluster import GeneCall
from gapit.errors import InputError
from gapit.typing_table import ScreenRow, parse_table

#: The table label for stdin input — the conventional ``-`` marker, also
#: how stdin renders in the JSON ``source`` array and error listings.
STDIN_LABEL = "-"


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
    except (OSError, UnicodeDecodeError) as exc:
        raise InputError(
            f"screen result is not valid UTF-8: {path}",
            code="SCREEN_TABLE_MALFORMED",
            context={"file": str(path)},
        ) from exc


def read_typing_input(paths: list[Path], stdin_text: str | None = None) -> TypingInput:
    """Parse and validate the designation input: every table's rows (the
    listed files plus, when ``stdin_text`` is given, the already-read stdin
    table) folded into per-FILE gene calls, with the single DATABASE they
    share."""
    rows: list[ScreenRow] = []
    labels: list[str] = []
    for path in paths:
        labels.append(str(path))
        rows.extend(parse_table(str(path), _read_table(path)))
    if stdin_text is not None:
        labels.append(STDIN_LABEL)
        rows.extend(parse_table(STDIN_LABEL, stdin_text))
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
