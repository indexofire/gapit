"""gapit.typing_result/1 — the `gapit typing` output document and renderers.

One result per (FILE, scheme): the designated phenotype with its
explainable breakdown. TSV/MD flatten each call to the seven reporting
columns (ambiguous calls render ``-`` and carry the candidate pair in
NOTES); JSON keeps the full :class:`gapit.typing_results.SchemeCall`
objects, the additive shape the inline engine produced before designation
moved to the command.
"""

from collections.abc import Iterable, Sequence
from datetime import datetime
from typing import Literal, NamedTuple

from pydantic import BaseModel, ConfigDict, Field

from gapit import __version__
from gapit.formats.json import ToolDocument, utc_timestamp
from gapit.formats.md import md_cell
from gapit.typing_results import PhenotypeScore, SchemeCall

COLUMNS = ("FILE", "SCHEME", "PHENOTYPE", "CONFIDENCE", "SCORE", "RUNNER_UP", "NOTES")


class TypingFileResult(NamedTuple):
    """One typed FILE and its scheme calls (scheme declaration order)."""

    file: str
    phenotypes: dict[str, SchemeCall]


class TypingFileDocument(BaseModel, frozen=True):
    """One typed FILE: scheme name → its call (the additive phenotypes
    object, kept from the inline-engine shape)."""

    file: str
    phenotypes: dict[str, SchemeCall]


class TypingResultDocument(BaseModel, frozen=True):
    """gapit.typing_result/1 — machine-readable designation output."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.typing_result/1"] = Field(
        default="gapit.typing_result/1", alias="schema"
    )
    tool: ToolDocument = ToolDocument()
    created_at: str
    source: list[str]
    db: str
    files: list[TypingFileDocument]


def _pair(entry: PhenotypeScore) -> str:
    return f"{entry.phenotype} ({entry.score:.4f})"


def _notes(call: SchemeCall) -> str:
    """The call's notes, ``;``-joined; an ambiguous call (null phenotype)
    appends its candidate pair so the flat table loses no signal."""
    notes = list(call.notes)
    if call.phenotype is None and call.ambiguous:
        notes.append("ambiguous: " + ", ".join(_pair(entry) for entry in call.ambiguous))
    return "; ".join(notes)


def typing_row(file: str, scheme: str, call: SchemeCall) -> tuple[str, ...]:
    """One result row's cells: the phenotype (``-`` when ambiguous, fallback
    strings verbatim — md.py's value semantics), the runner-up (``-`` when
    none, and also on an ambiguous call whose candidate pair NOTES
    carries), and the flattened notes."""
    called = call.runner_up if call.phenotype is not None else None
    return (
        file,
        scheme,
        call.phenotype or "-",
        call.confidence,
        f"{call.score:.4f}",
        "-" if called is None else _pair(called),
        _notes(call),
    )


def typing_tsv_preamble() -> str:
    """The header chunk of the typing TSV."""
    return "\t".join(COLUMNS) + "\n"


def typing_file_chunk(result: TypingFileResult) -> str:
    """One FILE's chunk: a row per scheme in declaration order — the
    streaming unit the use-case emits per typed file."""
    return "".join(
        "\t".join(typing_row(result.file, scheme, call)) + "\n"
        for scheme, call in result.phenotypes.items()
    )


def render_typing_tsv(results: Iterable[TypingFileResult]) -> str:
    """Render designation results as deterministic TSV."""
    return typing_tsv_preamble() + "".join(typing_file_chunk(result) for result in results)


def render_typing_result_json(
    source: Sequence[str], db: str, results: Iterable[TypingFileResult], *, now: datetime
) -> str:
    """Serialize designation results as gapit.typing_result/1 (indented,
    schema first)."""
    document = TypingResultDocument(
        created_at=utc_timestamp(now),
        source=list(source),
        db=db,
        files=[
            TypingFileDocument(file=result.file, phenotypes=result.phenotypes) for result in results
        ],
    )
    return document.model_dump_json(indent=2, by_alias=True)


def render_typing_md(
    source: Sequence[str], db: str, results: Iterable[TypingFileResult], *, now: datetime
) -> str:
    """Render designation results as deterministic Markdown: YAML
    frontmatter (tool, db, source tables) then one seven-column table."""
    files = list(results)
    lines: list[str] = [
        "---",
        "schema: gapit.typing_result/1",
        f"tool: gapit {__version__}",
        f"created_at: {utc_timestamp(now)}",
        f"db: {db}",
        "source:",
        *(f"  - {entry}" for entry in source),
        f"files: {len(files)}",
        "---",
        "",
        "# gapit typing results",
        "",
        "| " + " | ".join(COLUMNS) + " |",
        "|" + "---|" * len(COLUMNS),
    ]
    for result in files:
        for scheme, call in result.phenotypes.items():
            cells = typing_row(result.file, scheme, call)
            lines.append("| " + " | ".join(md_cell(cell) for cell in cells) + " |")
    return "\n".join(lines) + "\n"
