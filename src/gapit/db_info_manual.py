"""The man-page rendering of ``gapit db info`` (plain text, no ANSI).

A terminal reader wants the shape `man` taught everyone: a header rule with
the section title, uppercase section headings, an indented body, and hanging
indents for numbered references and bulleted caveats. The rendering is
deliberately **plain**: no colour, no box drawing, and a fixed width derived
from the terminal (falling back to 80) so a capture diffs cleanly and a pipe
still produces the byte-stable TSV instead (the ``db list`` precedent).

Only sections the document actually carries are printed — an uninstalled
provider has no PROVENANCE block, a user-built database without a sidecar has
no REFERENCES/CAVEATS/USAGE — so the page never lies with empty headings.
"""

import shutil
import textwrap

from gapit.db_info_ops import DbInfoDocument

_INDENT = " " * 7
_MIN_WIDTH = 60
_MAX_WIDTH = 100


def _width(width: int | None) -> int:
    if width is not None:
        return width
    columns = shutil.get_terminal_size((80, 24)).columns
    return max(_MIN_WIDTH, min(_MAX_WIDTH, columns))


def _heading(title: str) -> list[str]:
    """A blank line then the uppercase section title (lines joined by newline)."""
    return ["", title]


def _body(text: str, width: int, indent: str = _INDENT) -> str:
    """Wrap one paragraph into the man body indent."""
    return textwrap.fill(
        text,
        width=width,
        initial_indent=indent,
        subsequent_indent=indent,
        break_long_words=False,
        break_on_hyphens=False,
    )


def _labelled(label: str, value: str, width: int, pad: int = 14) -> str:
    """One ``Label  value`` row, wrapped with a hanging indent under the value."""
    column = _INDENT + " " * pad
    return textwrap.fill(
        value,
        width=width,
        initial_indent=f"{_INDENT}{label:<{pad - 1}} ",
        subsequent_indent=column,
        break_long_words=False,
        break_on_hyphens=False,
    )


def _header(name: str, width: int) -> str:
    """The classic three-column man header: name(1), title, name(1)."""
    page = f"{name.upper()}(1)"
    centre = "gapit database manual"
    gap = width - len(page) * 2 - len(centre)
    if gap < 2:  # too narrow for the centred title: keep the two page labels
        return f"{page}{' ' * max(2, width - len(page) * 2)}{page}"
    left = gap // 2
    return f"{page}{' ' * left}{centre}{' ' * (gap - left)}{page}"


def _synopsis(document: DbInfoDocument) -> tuple[str, ...]:
    """The commands this database actually participates in."""
    if document.kind == "cluster":
        return (
            f"gapit screen -d {document.name} [FILE]...",
            f"gapit db info -d {document.name}",
        )
    return (
        f"gapit screen -d {document.name} [FILE]... -o RESULT.tsv",
        "gapit typing RESULT.tsv",
        f"gapit db info -d {document.name}",
    )


def _status(document: DbInfoDocument) -> str:
    parts = [document.status]
    if document.records is not None:
        unit = "locus" if document.kind == "cluster" else "record"
        parts.append(f"{document.records} {unit}{'' if document.records == 1 else 's'}")
    parts.append(f"{document.kind} database")
    parts.append(document.dbtype)
    return " · ".join(parts)


def _provenance(document: DbInfoDocument, width: int) -> list[str]:
    rows: list[str] = []

    def add(label: str, value: str | None) -> None:
        if value:
            rows.append(_labelled(label, value, width))

    add("Vendor", document.vendor)
    add("Licence", document.license)
    add("Snapshot", document.snapshot)
    add("Typing", document.typing_schema)
    if document.build is not None:
        add("Built by", f"{document.build.tool} on {document.build.fetched_at}")
        tools: list[str] = []
        if document.build.makeblastdb_version:
            tools.append(
                document.build.makeblastdb_version.replace("makeblastdb: ", "makeblastdb ")
            )
        if document.build.minimap2_version:
            tools.append(f"minimap2 {document.build.minimap2_version}")
        add("Tooling", " · ".join(tools))
    add("Source", ", ".join(document.source_urls))
    if document.destination:
        add("Location", document.destination)
    return rows


def _references(document: DbInfoDocument, width: int) -> list[str]:
    """Numbered references with a hanging indent; doi/url/note stay attached."""
    rows: list[str] = []
    for index, reference in enumerate(document.references, start=1):
        marker = f"[{index}] "
        hanging = _INDENT + " " * len(marker)
        rows.append(
            textwrap.fill(
                reference.citation,
                width=width,
                initial_indent=_INDENT + marker,
                subsequent_indent=hanging,
                break_long_words=False,
                break_on_hyphens=False,
            )
        )
        for label, value in (("doi", reference.doi), ("url", reference.url)):
            if value:
                rows.append(f"{hanging}{label}: {value}")
        if reference.note:
            rows.append(
                textwrap.fill(
                    reference.note,
                    width=width,
                    initial_indent=hanging,
                    subsequent_indent=hanging,
                    break_long_words=False,
                    break_on_hyphens=False,
                )
            )
    return rows


def _usage(step: str, width: int) -> list[str]:
    """A usage step stays on one line; its ``# comment`` moves below it.

    Commands must never be wrapped mid-flag (a copied page has to stay
    runnable), so the aligned trailing comments the guidance carries become
    their own indented, wrapped lines.
    """
    command, _, comment = step.partition("#")
    lines = [f"{_INDENT}$ {command.rstrip()}"]
    note = comment.strip()
    if note:
        lines.append(
            textwrap.fill(
                f"# {note}",
                width=width,
                initial_indent=_INDENT + "  ",
                subsequent_indent=_INDENT + "  ",
                break_long_words=False,
                break_on_hyphens=False,
            )
        )
    return lines


def _bullets(items: tuple[str, ...], width: int) -> list[str]:
    return [
        textwrap.fill(
            item,
            width=width,
            initial_indent=f"{_INDENT}• ",
            subsequent_indent=_INDENT + "  ",
            break_long_words=False,
            break_on_hyphens=False,
        )
        for item in items
    ]


def render_db_manual(document: DbInfoDocument, width: int | None = None) -> str:
    """Render one ``gapit.dbinfo/1`` document as a plain man-style page."""
    columns = _width(width)
    lines: list[str] = [_header(document.name, columns)]

    lines += _heading("NAME")
    name_indent = _INDENT + " " * (len(document.name) + 3)
    lines.append(
        textwrap.fill(
            f"— {document.description}",
            width=columns,
            initial_indent=_INDENT + document.name + " ",
            subsequent_indent=name_indent,
            break_long_words=False,
            break_on_hyphens=False,
        )
    )

    lines += _heading("SYNOPSIS")
    lines += [_body(f"$ {command}", columns) for command in _synopsis(document)]

    if document.description or document.note:
        lines += _heading("DESCRIPTION")
        if document.note:
            lines.append(_body(document.note, columns))

    lines += _heading("STATUS")
    lines.append(_body(_status(document), columns))

    provenance = _provenance(document, columns)
    if provenance:
        lines += _heading("PROVENANCE")
        lines += provenance

    if document.references:
        lines += _heading("REFERENCES")
        lines += _references(document, columns)

    if document.caveats:
        lines += _heading("CAVEATS")
        lines += _bullets(document.caveats, columns)

    if document.usage:
        lines += _heading("USAGE")
        for step in document.usage:
            lines += _usage(step, columns)

    lines += _heading("SEE ALSO")
    lines.append(
        _body(
            "gapit db list(1), gapit db search(1), gapit screen(1), gapit typing(1), "
            "docs/databases.md",
            columns,
        )
    )
    return "\n".join(lines).rstrip() + "\n"
