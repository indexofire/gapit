"""The `gapit summary` command: matrix over abricate-format report tables.

Lives outside cli.py to keep that module small; cli.py registers it via
``register_summary_command``.
"""

from datetime import UTC, datetime
from pathlib import Path
from typing import Annotated

import typer

from gapit.dispatch import dispatch
from gapit.formats.summary import format_summary_tsv, render_summary_json, render_summary_md
from gapit.screening import OutputFormat, usage_fail
from gapit.summary import SummaryParams, build_summary


def summary_command(
    files: Annotated[
        list[Path] | None,
        typer.Argument(help="Abricate-format report file(s) to summarize."),
    ] = None,
    identity: Annotated[
        bool,
        typer.Option("--identity", help="Cells show %IDENTITY instead of %COVERAGE."),
    ] = False,
    nopath: Annotated[
        bool,
        typer.Option("--nopath", help="Basename row keys (FILE values / input filenames)."),
    ] = False,
    quiet: Annotated[bool, typer.Option("--quiet", help="Silence stderr diagnostics.")] = False,
    output_format: Annotated[
        OutputFormat, typer.Option("--format", help="Output format.")
    ] = OutputFormat.tsv,
) -> None:
    """Summarize report table(s) into a gene presence/absence matrix."""

    def run() -> None:
        if not files:
            usage_fail("summary needs >= 1 report file(s)")
        if identity and not quiet:
            typer.echo("Using %IDENTITY for the summary table instead of %COVERAGE", err=True)
        params = SummaryParams(identity=identity, nopath=nopath)

        def warn(message: str) -> None:
            if not quiet:
                typer.echo(f"WARNING: {message}", err=True)

        matrix = build_summary(list(files), params, warn=warn)
        now = datetime.now(UTC)
        if output_format is OutputFormat.json:
            typer.echo(render_summary_json(matrix, now=now), nl=False)
        elif output_format is OutputFormat.md:
            typer.echo(render_summary_md(matrix, now=now), nl=False)
        else:
            typer.echo(format_summary_tsv(matrix, csv=output_format is OutputFormat.csv), nl=False)

    dispatch(run)


def register_summary_command(app: typer.Typer) -> None:
    """Attach the summary command to the CLI app."""
    app.command("summary")(summary_command)
