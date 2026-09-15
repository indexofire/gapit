"""gaita command-line interface (typer entrypoint)."""

from typing import Annotated

import typer

from gaita import __version__

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
