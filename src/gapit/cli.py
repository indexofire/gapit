"""gapit command-line interface (typer entrypoint)."""

import json
from collections.abc import Callable
from pathlib import Path
from typing import Annotated

import typer
from pydantic import BaseModel

from gapit import __version__, config, db
from gapit.cmd_db import register_db_command
from gapit.cmd_screen import register_screen_command
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
from gapit.screening import usage_fail

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


register_screen_command(app)
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
