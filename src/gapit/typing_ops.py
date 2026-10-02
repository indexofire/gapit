"""The `gapit typing` use-case: designation from screen result tables.

Stage 2 of the two-stage pipeline (rightsholder design): ``gapit screen
-o result.tsv --db NAME`` detects genes; this use-case reads the table(s),
resolves NAME from the DATABASE column under --datadir (the default
resolution — bundled databases still materialize on first use), loads the
database's typing.json, and runs the shared gene-path engine
(:func:`gapit.typing_gene.evaluate_gene_calls`) over each FILE's folded
gene calls. Cluster databases refuse with a typed error: their typing is
integrated into ``gapit screen`` (gapit.cluster/1 emits it directly).
"""

from datetime import UTC, datetime
from pathlib import Path

import typer

from gapit.bundled import resolve_screen_datadir
from gapit.cluster import load_typing
from gapit.engines import Emit, OutputFormat
from gapit.errors import DatabaseError
from gapit.formats.typing_result import (
    TypingFileResult,
    render_typing_md,
    render_typing_result_json,
    typing_file_chunk,
    typing_tsv_preamble,
)
from gapit.screening import find_database
from gapit.typing_gene import evaluate_gene_calls
from gapit.typing_input import STDIN_LABEL, read_typing_input


def _resolve_document(table_database: str, datadir: Path | None, quiet: bool):
    """The typing document of the database every row screened against —
    refusing cluster databases (integrated typing) and untyped ones."""
    database = find_database(
        resolve_screen_datadir(datadir, table_database), table_database, quiet=quiet
    )
    if database.kind == "cluster":
        raise DatabaseError(
            f"database {table_database!r} is a cluster database; cluster typing is integrated"
            f" into screening — `gapit screen --db {table_database}` already emits the phenotype",
            code="TYPING_CLUSTER_DB",
            context={"db": table_database},
        )
    document = load_typing(database)
    if document is None:
        raise DatabaseError(
            f"database {table_database!r} carries no typing scheme"
            " (only databases built with `db build --typing` designate)",
            code="TYPING_NO_SCHEME",
            context={"db": table_database},
        )
    return document


def run_typing(
    paths: list[Path],
    datadir: Path | None,
    output_format: OutputFormat,
    quiet: bool,
    emit: Emit | None = None,
    stdin_text: str | None = None,
) -> str:
    """Type every FILE of the screen result table(s) and render; returns the
    full output for the caller to echo and, when ``emit`` is given, hands
    each rendered chunk to it as soon as it exists (tsv streams per typed
    file; md/json are single documents). The table(s) arrive as ``paths``
    and/or as ``stdin_text`` — the already-read piped table, rendered under
    the ``-`` source label."""
    source = [str(path) for path in paths]
    if stdin_text is not None:
        source.append(STDIN_LABEL)
    table = read_typing_input(paths, stdin_text)
    document = _resolve_document(table.database, datadir, quiet)
    if not quiet:
        typer.echo(f"Typing {len(table.files)} file(s) against {table.database}", err=True)
    results = [
        TypingFileResult(file=entry.file, phenotypes=evaluate_gene_calls(entry.calls, document))
        for entry in table.files
    ]
    now = datetime.now(UTC)
    chunks: list[str] = []

    def sink(chunk: str) -> None:
        if chunk:
            chunks.append(chunk)
            if emit is not None:
                emit(chunk)

    if output_format is OutputFormat.json:
        sink(render_typing_result_json(source, table.database, results, now=now))
    elif output_format is OutputFormat.md:
        sink(render_typing_md(source, table.database, results, now=now))
    else:
        sink(typing_tsv_preamble())
        for result in results:
            sink(typing_file_chunk(result))
    return "".join(chunks)
