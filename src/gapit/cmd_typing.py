"""The `gapit typing` command: designation from screen result tables.

Lives outside cli.py to keep that module small; cli.py registers it via
``register_typing_command``. The second stage of the designation pipeline:
``gapit screen -o result.tsv --db NAME`` detects genes, then
``gapit typing result.tsv`` resolves NAME from the table's DATABASE column
and designates (the use-case lives in gapit.typing_ops). The table also
arrives on stdin: ``gapit screen 1.fna --db NAME | gapit typing`` types the
piped output, and the explicit ``-`` marker reads stdin even under a
terminal (until EOF).
"""

import sys
from pathlib import Path

import click

from gapit.dispatch import EnumChoice, datadir_option, dispatch, output_target, quiet_option
from gapit.engines import OutputFormat
from gapit.errors import InputError, usage_fail
from gapit.typing_input import STDIN_LABEL
from gapit.typing_ops import run_typing


def _stdin_is_tty() -> bool:
    """Whether stdin is a terminal — the injectable probe behind the
    bare-help guard (a CliRunner stdin is never a tty)."""
    return sys.stdin.isatty()


def _read_stdin() -> str:
    """The piped screen table's text; non-UTF-8 bytes are an InputError."""
    try:
        return sys.stdin.read()
    except UnicodeDecodeError as exc:
        raise InputError(
            "screen result piped to stdin is not valid UTF-8",
            code="SCREEN_TABLE_MALFORMED",
            context={"file": STDIN_LABEL},
        ) from exc


def _stdin_text(files: list[Path]) -> str | None:
    """Stdin's table text when stdin is the input: the explicit ``-``
    marker (alone — mixing it with file arguments is a usage error) or the
    argument-less piped invocation."""
    if any(str(path) == STDIN_LABEL for path in files):
        if len(files) > 1:
            usage_fail("type either '-' (stdin) or screen result file(s), not both")
        return _read_stdin()
    if files:
        return None
    return _read_stdin()


@click.argument(
    "files",
    nargs=-1,
    required=False,
    type=click.Path(path_type=Path),
    help="gapit/abricate screen result TSV/CSV file(s) to type ('-' = stdin).",
)
@datadir_option()
@click.option(
    "--format",
    "-f",
    "output_format",
    type=EnumChoice(OutputFormat),
    default=OutputFormat.tsv,
    help="Output format.",
)
@quiet_option()
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    help=(
        "Write the typing results to PATH instead of stdout (truncates; stdout"
        " then carries no data)."
    ),
)
@click.pass_context
def typing_command(
    ctx: click.Context,
    files: tuple[Path, ...],
    datadir: Path | None,
    output_format: OutputFormat,
    quiet: bool,
    output: Path | None,
) -> None:
    """Designate phenotypes from screen result table(s) (typed gene databases)."""
    args = list(files or [])
    # Bare invocation at a terminal keeps the help-on-bare policy; with a
    # pipe attached the table arrives on stdin. Runs before dispatch so the
    # help path exits cleanly — never through the error-envelope handler.
    if not args and _stdin_is_tty():
        click.echo(ctx.get_help())
        ctx.exit(2)

    def run() -> None:
        if output_format is OutputFormat.csv:
            usage_fail("--format csv is not available for typing (use tsv, json, or md)")
        stdin_text = _stdin_text(args)
        paths = args if stdin_text is None else []
        with output_target(output) as deliver:
            # run_typing emits its output through `deliver`; echoing the
            # return value here would duplicate every byte.
            run_typing(paths, datadir, output_format, quiet, emit=deliver, stdin_text=stdin_text)

    dispatch(run)


def register_typing_command(app: click.Group) -> None:
    """Attach the typing command to the CLI app (bare invocation at a
    terminal prints help; piped stdin types without arguments)."""
    app.command("typing")(typing_command)
