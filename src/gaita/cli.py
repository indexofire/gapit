"""gaita command-line interface (typer entrypoint)."""

import json
from collections.abc import Callable
from pathlib import Path
from typing import Annotated

import typer

from gaita import __version__, config, db
from gaita.errors import GaitaError

app = typer.Typer(
    name="gaita",
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
        typer.echo(f"gaita {__version__}")
        raise typer.Exit()
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise typer.Exit()


Datadir = Annotated[
    Path | None,
    typer.Option(
        "--datadir",
        help="Database directory (default: $GAITA_DATADIR, then ~/.local/share/gaita/db).",
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
            "schema": "gaita.list/1",
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
