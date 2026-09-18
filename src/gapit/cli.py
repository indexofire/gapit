"""gapit command-line interface (typer entrypoint)."""

import json
from collections.abc import Callable
from pathlib import Path
from typing import Annotated

import typer
from pydantic import BaseModel

from gapit import __version__, config, db
from gapit.cmd_db import register_db_command
from gapit.cmd_summary import register_summary_command
from gapit.errors import ErrorEnvelope, GapitError, render_error
from gapit.formats.json import (
    ListDocument,
    ListEntryDocument,
    ReadsDocument,
    ReportDocument,
    VersionDocument,
)
from gapit.formats.summary import SummaryDocument
from gapit.mcp import register_mcp_command
from gapit.reads import ReadTypeEnum
from gapit.screening import OutputFormat, run_screen, run_screen_reads, usage_fail

app = typer.Typer(
    name="gapit",
    help="Mass screening of contigs for antimicrobial resistance and virulence genes.",
    no_args_is_help=True,
    add_completion=True,
)


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    show_version: Annotated[
        bool | None,
        typer.Option("--version", help="Show version and exit."),
    ] = None,
    as_json: Annotated[
        bool,
        typer.Option("--json", help="With --version: emit gapit.version/1 JSON."),
    ] = False,
) -> None:
    """Mass screening of contigs for AMR and virulence genes."""
    if show_version:
        if as_json:
            typer.echo(VersionDocument(version=__version__).model_dump_json(by_alias=True))
        else:
            typer.echo(f"gapit {__version__}")
        raise typer.Exit()
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise typer.Exit()


Datadir = Annotated[
    Path | None,
    typer.Option(
        "--datadir",
        help="Database directory (default: $GAPIT_DATADIR, then ~/.local/share/gapit/db).",
    ),
]


def _dispatch(action: Callable[[], None]) -> None:
    """Run a command body; any failure renders the gapit.error/1 envelope on
    stderr and exits with the documented code (UNEXPECTED/1 for non-GapitError)."""
    try:
        action()
    except Exception as exc:
        typer.echo(render_error(exc), err=True)
        raise typer.Exit(code=exc.exit_code if isinstance(exc, GapitError) else 1) from exc


def _list(datadir: Path | None, as_json: bool) -> None:
    infos = db.list_databases(config.resolve_datadir(datadir), setupdb=False)
    if as_json:
        document = ListDocument(
            databases=[
                ListEntryDocument(
                    name=info.name,
                    sequences=info.n_sequences,
                    dbtype=info.dbtype,
                    date=info.date,
                )
                for info in infos
            ]
        )
        typer.echo(document.model_dump_json(indent=2, by_alias=True))
        return
    typer.echo("DATABASE\tSEQUENCES\tDBTYPE\tDATE")
    for info in infos:
        typer.echo(f"{info.name}\t{info.n_sequences}\t{info.dbtype}\t{info.date}")


def _setupdb(datadir: Path | None, debug: bool) -> None:
    infos = db.list_databases(config.resolve_datadir(datadir), setupdb=True, debug=debug)
    for info in infos:
        typer.echo(
            f"Indexed {info.name} ({info.n_sequences} sequences, {info.dbtype})",
            err=True,
        )


@app.command("list")
def list_dbs(
    datadir: Datadir = None,
    as_json: Annotated[
        bool,
        typer.Option("--json", help="Print machine-readable JSON instead of a table."),
    ] = False,
) -> None:
    """List installed databases (abricate --list compatible)."""
    _dispatch(lambda: _list(datadir, as_json))


@app.command("setupdb")
def setupdb(
    datadir: Datadir = None,
    debug: Annotated[
        bool,
        typer.Option("--debug", help="Echo external command lines to stderr."),
    ] = False,
) -> None:
    """Build BLAST indices for all databases under the datadir."""
    _dispatch(lambda: _setupdb(datadir, debug))


@app.command("screen")
def screen(
    files: Annotated[
        list[Path] | None,
        typer.Argument(help="Input FASTA/GBK/EMBL contig file(s) to screen."),
    ] = None,
    r1: Annotated[
        str | None,
        typer.Option("--r1", help="Comma-separated FASTQ R1 file(s), one per lane (reads mode)."),
    ] = None,
    r2: Annotated[
        str | None,
        typer.Option("--r2", help="Comma-separated mate FASTQ file(s); must match --r1 count."),
    ] = None,
    read_type: Annotated[
        ReadTypeEnum,
        typer.Option(
            "--read-type",
            help="minimap2 preset for reads mode (sr, map-ont, map-hifi).",
        ),
    ] = ReadTypeEnum.sr,
    min_breadth: Annotated[
        float,
        typer.Option("--min-breadth", help="Reads mode: minimum %breadth for presence."),
    ] = 90.0,
    db: Annotated[
        str, typer.Option("--db", help="Database to screen against (datadir subdir).")
    ] = "ncbi",
    datadir: Datadir = None,
    minid: Annotated[
        float, typer.Option("--minid", help="Minimum %identity, 0 < x <= 100.")
    ] = 80.0,
    mincov: Annotated[
        float, typer.Option("--mincov", help="Minimum %coverage, 0 <= x <= 100.")
    ] = 80.0,
    threads: Annotated[int, typer.Option("--threads", help="BLAST worker threads.")] = 1,
    jobs: Annotated[
        int,
        typer.Option(
            "--jobs",
            help=(
                "Screen N input files concurrently (gapit extension; output order is"
                " always input order). Each worker runs its own BLAST against the"
                " shared db index, which BLAST mmaps — concurrent readers are fine."
            ),
        ),
    ] = 1,
    fofn: Annotated[
        Path | None,
        typer.Option("--fofn", help="File of filenames; replaces the positional FILEs."),
    ] = None,
    quiet: Annotated[bool, typer.Option("--quiet", help="Silence stderr diagnostics.")] = False,
    csv_flag: Annotated[bool, typer.Option("--csv", help="Compat alias for --format csv.")] = False,
    noheader: Annotated[bool, typer.Option("--noheader", help="Suppress the header row.")] = False,
    nopath: Annotated[bool, typer.Option("--nopath", help="Basename the FILE column.")] = False,
    debug: Annotated[bool, typer.Option("--debug", help="Verbose stderr diagnostics.")] = False,
    output_format: Annotated[
        OutputFormat | None,
        typer.Option("--format", help="Output format (reads mode defaults to json)."),
    ] = None,
) -> None:
    """Screen contig files or FASTQ reads (R1 and R2 comma-lists, one lane
    each) for AMR/virulence genes."""

    def run() -> None:
        if (r1 is not None or r2 is not None) and files:
            usage_fail("--r1/--r2 and positional contig FILEs are mutually exclusive")
        if r2 is not None and r1 is None:
            usage_fail("--r2 requires --r1")
        if r1 is not None or r2 is not None:
            run_screen_reads(
                r1 or "",
                r2,
                db,
                datadir,
                read_type,
                min_breadth,
                threads,
                output_format,
                quiet,
                debug,
            )
        else:
            run_screen(
                files,
                db,
                datadir,
                minid,
                mincov,
                threads,
                jobs,
                fofn,
                quiet,
                csv_flag,
                noheader,
                nopath,
                debug,
                output_format or OutputFormat.tsv,
            )

    _dispatch(run)


register_summary_command(app)
register_db_command(app)
register_mcp_command(app)


_SCHEMA_MODELS: dict[str, type[BaseModel]] = {
    "report": ReportDocument,
    "reads": ReadsDocument,
    "summary": SummaryDocument,
    "list": ListDocument,
    "error": ErrorEnvelope,
    "version": VersionDocument,
}


@app.command("schema")
def schema(
    name: Annotated[
        str,
        typer.Argument(
            help="Document to introspect: report, reads, summary, list, error, or version."
        ),
    ],
) -> None:
    """Print the JSON Schema of a gapit output document."""

    def run() -> None:
        model = _SCHEMA_MODELS.get(name)
        if model is None:
            usage_fail(f"unknown schema name: {name} (choose from: {', '.join(_SCHEMA_MODELS)})")
        typer.echo(json.dumps(model.model_json_schema(by_alias=True), indent=2))

    _dispatch(run)
