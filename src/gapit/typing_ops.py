"""The `gspec typing` use-case: designation from screen result tables.

Stage 2 of the two-stage pipeline (rightsholder design): ``gspec screen
-o result.tsv --db NAME`` detects genes (contig mode, or reads mode's
``#SAMPLE``/``BREADTH%`` table); this use-case reads the table(s), resolves
NAME from the DATABASE column under --datadir (the default resolution —
bundled databases still materialize on first use), loads the database's
typing.json, and runs the shared gene-path engine
(:func:`gapit.typing_gene.evaluate_gene_calls`) over each FILE's folded
gene calls. Cluster databases refuse with a typed error: their typing is
integrated into ``gspec screen``.  Cluster JSON results (``gapit.cluster/2``
or ``gapit.locus_map/1``) are accepted as a concise-formatting path: their
already-evaluated ``phenotypes`` map is rendered as the same
``gapit.typing_result/1`` document, without re-running the aligner.
"""

from datetime import UTC, datetime
from pathlib import Path

import click

from gapit.bundled import resolve_screen_datadir
from gapit.cluster import load_typing
from gapit.db import Database
from gapit.dbcodec import decode_seqid
from gapit.engines import Emit, OutputFormat
from gapit.errors import DatabaseError
from gapit.fasta import iter_fasta
from gapit.formats.typing_result import (
    TypingFileResult,
    render_typing_md,
    render_typing_result_json,
    typing_file_chunk,
    typing_tsv_preamble,
)
from gapit.screening import find_database
from gapit.serotype_engine import evaluate_serotype
from gapit.serotype_models import load_serotype, validate_serotype_references
from gapit.typing_cluster_input import read_cluster_typing_input
from gapit.typing_gene import evaluate_gene_calls
from gapit.typing_input import STDIN_LABEL, read_typing_input


def _resolve_database(table_database: str, datadir: Path | None, quiet: bool) -> Database:
    """The database every row screened against, refusing cluster databases."""
    database = find_database(
        resolve_screen_datadir(datadir, table_database), table_database, quiet=quiet
    )
    if database.kind == "cluster":
        raise DatabaseError(
            f"database {table_database!r} is a cluster database; cluster typing is integrated"
            f" into screening — `gspec screen --db {table_database}` already emits the phenotype",
            code="TYPING_CLUSTER_DB",
            context={"db": table_database},
        )
    return database


def _typing_document(database: Database):
    document = load_typing(database)
    if document is None:
        raise DatabaseError(
            f"database {database.name!r} carries no typing scheme"
            " (only databases built with `db build --typing` designate)",
            code="TYPING_NO_SCHEME",
            context={"db": database.name},
        )
    return document


def _database_gene_ids(database: Database) -> frozenset[str]:
    return frozenset(
        decode_seqid(record.id, database.name).gene
        for record in iter_fasta(database.sequences_path)
    )


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

    def render(db_name: str, results: list[TypingFileResult], *, now: datetime) -> str:
        chunks: list[str] = []

        def sink(chunk: str) -> None:
            if chunk:
                chunks.append(chunk)
                if emit is not None:
                    emit(chunk)

        if output_format is OutputFormat.json:
            sink(render_typing_result_json(source, db_name, results, now=now))
        elif output_format is OutputFormat.md:
            sink(render_typing_md(source, db_name, results, now=now))
        else:
            sink(typing_tsv_preamble())
            for result in results:
                sink(typing_file_chunk(result))
        return "".join(chunks)

    cluster_input = read_cluster_typing_input(paths, stdin_text)
    if cluster_input is not None:
        database, results = cluster_input
        if not quiet:
            click.echo(
                f"Typing {len(results)} file(s) against {database} (cluster result)",
                err=True,
            )
        return render(database, results, now=datetime.now(UTC))

    table = read_typing_input(paths, stdin_text)
    database = _resolve_database(table.database, datadir, quiet)
    serotype_document = load_serotype(database.path)
    if serotype_document is not None:
        validate_serotype_references(serotype_document, _database_gene_ids(database))
        if not quiet:
            click.echo(
                f"Typing {len(table.files)} file(s) against {table.database} (serotype)",
                err=True,
            )
        results = [
            TypingFileResult(
                file=entry.file,
                genes=tuple(sorted(entry.calls)),
                phenotypes={"serotype": evaluate_serotype(entry.calls, serotype_document)},
            )
            for entry in table.files
        ]
    else:
        document = _typing_document(database)
        if not quiet:
            click.echo(f"Typing {len(table.files)} file(s) against {table.database}", err=True)
        results = [
            TypingFileResult(
                file=entry.file,
                genes=tuple(sorted(entry.calls)),
                phenotypes=evaluate_gene_calls(entry.calls, document),
            )
            for entry in table.files
        ]
    return render(table.database, results, now=datetime.now(UTC))
