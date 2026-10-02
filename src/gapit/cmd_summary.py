"""The `gapit summary` command: matrix over abricate-format report tables.

Lives outside cli.py to keep that module small; cli.py registers it via
``register_summary_command``. The report table also arrives on stdin:
``gapit screen -d ecoli_dec *.fna | gapit summary`` summarizes the piped
output (a batched multi-FILE table is one input, so dutch mode applies),
and the explicit ``-`` marker reads stdin even under a terminal (until
EOF). A bare invocation at a terminal keeps the help-on-bare policy.
"""

import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Annotated

import typer

from gapit.dispatch import dispatch
from gapit.errors import InputError, usage_fail
from gapit.formats.summary import format_summary_tsv, render_summary_json, render_summary_md
from gapit.screening import OutputFormat
from gapit.summary import STDIN_LABEL, SummaryParams, build_summary


def _stdin_is_tty() -> bool:
    """Whether stdin is a terminal — the injectable probe behind the
    bare-help guard (a CliRunner stdin is never a tty). Duplicated from
    cmd_typing so each module attribute stays its tests' monkeypatch seam."""
    return sys.stdin.isatty()


def _read_stdin() -> str:
    """The piped report table's text; non-UTF-8 bytes are an InputError."""
    try:
        return sys.stdin.read()
    except UnicodeDecodeError as exc:
        raise InputError(
            "report piped to stdin is not valid UTF-8",
            code="SUMMARY_MALFORMED",
            context={"file": STDIN_LABEL},
        ) from exc


def _stdin_text(files: list[Path]) -> str | None:
    """Stdin's table text when stdin is the input: the explicit ``-``
    marker (alone — mixing it with file arguments is a usage error) or the
    argument-less piped invocation."""
    if any(str(path) == STDIN_LABEL for path in files):
        if len(files) > 1:
            usage_fail("summarize either '-' (stdin) or report file(s), not both")
        return _read_stdin()
    if files:
        return None
    return _read_stdin()


def summary_command(
    ctx: typer.Context,
    files: Annotated[
        list[Path] | None,
        typer.Argument(help="Abricate-format report file(s) to summarize ('-' = stdin)."),
    ] = None,
    identity: Annotated[
        bool,
        typer.Option(
            "--identity", "-i", help="Cells show %IDENTITY values (default: +/- presence)."
        ),
    ] = False,
    coverage: Annotated[
        bool,
        typer.Option(
            "--coverage", "-c", help="Cells show %COVERAGE values (default: +/- presence)."
        ),
    ] = False,
    nopath: Annotated[
        bool,
        typer.Option("--nopath", "-p", help="Basename row keys (FILE values / input filenames)."),
    ] = False,
    quiet: Annotated[
        bool, typer.Option("--quiet", "-q", help="Silence stderr diagnostics.")
    ] = False,
    output_format: Annotated[
        OutputFormat, typer.Option("--format", "-f", help="Output format.")
    ] = OutputFormat.tsv,
) -> None:
    """Summarize report table(s) into a gene presence/absence matrix."""
    args = list(files or [])
    # Bare invocation at a terminal keeps the help-on-bare policy; with a
    # pipe attached the table arrives on stdin. Runs before dispatch so the
    # help path exits cleanly — never through the error-envelope handler.
    if not args and _stdin_is_tty():
        typer.echo(ctx.get_help())
        raise typer.Exit(code=2)

    def run() -> None:
        stdin_text = _stdin_text(args)
        paths = args if stdin_text is None else []
        if identity and coverage and not quiet:
            typer.echo("Using %IDENTITY/%COVERAGE per hit for the summary table", err=True)
        elif identity and not quiet:
            typer.echo("Using %IDENTITY for the summary table instead of +/- presence", err=True)
        elif coverage and not quiet:
            typer.echo("Using %COVERAGE for the summary table instead of +/- presence", err=True)
        params = SummaryParams(identity=identity, coverage=coverage, nopath=nopath)

        def warn(message: str) -> None:
            if not quiet:
                typer.echo(f"WARNING: {message}", err=True)

        matrix = build_summary(paths, params, warn=warn, stdin_text=stdin_text)
        now = datetime.now(UTC)
        if output_format is OutputFormat.json:
            typer.echo(render_summary_json(matrix, now=now), nl=False)
        elif output_format is OutputFormat.md:
            typer.echo(render_summary_md(matrix, now=now), nl=False)
        else:
            typer.echo(format_summary_tsv(matrix, csv=output_format is OutputFormat.csv), nl=False)

    dispatch(run)


def register_summary_command(app: typer.Typer) -> None:
    """Attach the summary command to the CLI app (bare invocation at a
    terminal prints help; piped stdin summarizes without arguments)."""
    app.command("summary")(summary_command)
