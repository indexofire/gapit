"""gapit command-line interface (typer entrypoint)."""

import enum
import json
from collections.abc import Callable
from pathlib import Path
from typing import Annotated, NoReturn

import typer

from gapit import __version__, config, db
from gapit.blast import screen_file
from gapit.errors import DatabaseError, GaitaError, InputError
from gapit.formats.tsv import format_tsv
from gapit.report import Report, ScreeningParams

app = typer.Typer(
    name="gapit",
    help="Mass screening of contigs for antimicrobial resistance and virulence genes.",
    no_args_is_help=True,
    add_completion=False,
)


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    show_version: Annotated[
        bool | None,
        typer.Option("--version", help="Show version and exit."),
    ] = None,
) -> None:
    """Mass screening of contigs for AMR and virulence genes."""
    if show_version:
        typer.echo(f"gapit {__version__}")
        raise typer.Exit()
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise typer.Exit()


Datadir = Annotated[
    Path | None,
    typer.Option(
        "--datadir",
        help="Database directory (default: $GAITA_DATADIR, then ~/.local/share/gapit/db).",
    ),
]


def _dispatch(action: Callable[[], None]) -> None:
    """Run a command body, mapping GaitaError to `ERROR: ...` on stderr + its exit code."""
    try:
        action()
    except GaitaError as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=exc.exit_code) from exc


def _list(datadir: Path | None, as_json: bool) -> None:
    infos = db.list_databases(config.resolve_datadir(datadir), setupdb=False)
    if as_json:
        payload = {
            "schema": "gapit.list/1",
            "databases": [
                {
                    "name": info.name,
                    "sequences": info.n_sequences,
                    "dbtype": info.dbtype,
                    "date": info.date,
                }
                for info in infos
            ],
        }
        typer.echo(json.dumps(payload, indent=2))
        return
    typer.echo("DATABASE\tSEQUENCES\tDBTYPE\tDATE")
    for info in infos:
        typer.echo(f"{info.name}\t{info.n_sequences}\t{info.dbtype}\t{info.date}")


def _setupdb(datadir: Path | None) -> None:
    infos = db.list_databases(config.resolve_datadir(datadir), setupdb=True)
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
def setupdb(datadir: Datadir = None) -> None:
    """Build BLAST indices for all databases under the datadir."""
    _dispatch(lambda: _setupdb(datadir))


class OutputFormat(enum.Enum):
    """Screen output formats (JSON/Markdown arrive in Phase 4)."""

    tsv = "tsv"
    csv = "csv"


def _usage_fail(message: str) -> NoReturn:
    """Report a usage error on stderr and exit 2."""
    typer.echo(f"ERROR: {message}", err=True)
    raise typer.Exit(code=2)


def _resolve_inputs(files: list[Path] | None, fofn: Path | None) -> list[Path]:
    """Input files: --fofn (lines stripped, empties dropped) REPLACES positionals."""
    if fofn is not None:
        if not fofn.is_file():
            raise InputError(f"--fofn file not found: {fofn}", code="INPUT_NOT_FOUND")
        inputs = [
            Path(line.strip())
            for line in fofn.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    elif files:
        inputs = list(files)
    else:
        _usage_fail("no input files given (positional FILEs or --fofn)")
    for path in inputs:
        if not path.is_file():
            raise InputError(f"input file not found or unreadable: {path}", code="INPUT_NOT_FOUND")
    return inputs


def _find_database(datadir: Path, name: str) -> db.Database:
    """Look up a database by name under the datadir; unknown names list what exists."""
    databases = db.discover_databases(datadir)
    for database in databases:
        if database.name == name:
            return database
    available = ", ".join(entry.name for entry in databases) or "(none)"
    raise DatabaseError(
        f"Database {name} is not in {datadir}. Available: {available}",
        code="DATABASE_NOT_FOUND",
    )


def _screen(
    files: list[Path] | None,
    db_name: str,
    datadir: Path | None,
    minid: float,
    mincov: float,
    threads: int,
    fofn: Path | None,
    quiet: bool,
    csv_flag: bool,
    noheader: bool,
    nopath: bool,
    debug: bool,
    output_format: OutputFormat,
) -> None:
    """Screen each input file in order; buffer reports; print once at the end."""
    if not 0.0 < minid <= 100.0:
        _usage_fail(f"--minid must be in (0, 100]: got {minid}")
    if not 0.0 <= mincov <= 100.0:
        _usage_fail(f"--mincov must be in [0, 100]: got {mincov}")
    if threads < 1:
        _usage_fail(f"--threads must be >= 1: got {threads}")
    inputs = _resolve_inputs(files, fofn)
    params = ScreeningParams(db=db_name, minid=minid, mincov=mincov, threads=threads)
    database = _find_database(config.resolve_datadir(datadir), db_name)
    reports: list[Report] = []
    for path in inputs:
        if not quiet:
            typer.echo(f"Processing: {path}", err=True)
        report = screen_file(path, database, params, debug=debug)
        if not quiet:
            typer.echo(f"Found {len(report.hits)} genes in {path}", err=True)
        reports.append(report)
    as_csv = csv_flag or output_format is OutputFormat.csv
    typer.echo(format_tsv(reports, csv=as_csv, noheader=noheader, nopath=nopath), nl=False)


@app.command("screen")
def screen(
    files: Annotated[
        list[Path] | None,
        typer.Argument(help="Input FASTA/GBK/EMBL file(s) to screen."),
    ] = None,
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
        OutputFormat, typer.Option("--format", help="Output format.")
    ] = OutputFormat.tsv,
) -> None:
    """Screen contig files for AMR/virulence genes (abricate-compatible TSV)."""
    _dispatch(
        lambda: _screen(
            files,
            db,
            datadir,
            minid,
            mincov,
            threads,
            fofn,
            quiet,
            csv_flag,
            noheader,
            nopath,
            debug,
            output_format,
        )
    )
