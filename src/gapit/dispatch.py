"""Shared CLI plumbing: error-envelope dispatch, the ``--datadir`` option,
and the ``--output`` delivery target (stdout chunks vs a lazily opened
file).

Imported by every command module (cli.py, cmd_*.py); imports nothing from
gapit except :mod:`gapit.errors`, so it can never participate in a cycle.
"""

from collections.abc import Callable, Generator
from contextlib import contextmanager
from pathlib import Path
from typing import IO, Annotated, TextIO

import typer

from gapit.errors import GapitError, render_error

Datadir = Annotated[
    Path | None,
    typer.Option(
        "--datadir",
        "-D",
        help="Database directory (default: $GAPIT_DATADIR, then ~/.local/share/gapit/db).",
    ),
]


class OutputFileWriter:
    """Chunk sink behind ``screen --output``: the file opens (truncating)
    on the first chunk — a run that fails before rendering anything leaves
    no empty artifact — and every chunk is flushed, so an in-batch failure
    still persists the already-streamed prefix. The caller closes it."""

    def __init__(self, path: Path) -> None:
        self._path = path
        self._handle: TextIO | None = None

    def __call__(self, chunk: str) -> None:
        handle: IO[str] = self._handle or self._path.open("w", encoding="utf-8", newline="\n")
        self._handle = handle
        handle.write(chunk)
        handle.flush()

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()


@contextmanager
def output_target(path: Path | None) -> Generator[Callable[[str], None]]:
    """Where rendered output goes: stdout (echoed per chunk, no trailing
    newline) or the file behind ``--output`` (stdout then carries no data).
    Single-document engines call the delivered callable once; streaming
    engines pass it as the use-case ``emit`` sink."""
    if path is None:
        yield lambda chunk: typer.echo(chunk, nl=False)
        return
    writer = OutputFileWriter(path)
    try:
        yield writer
    finally:
        writer.close()


def dispatch(action: Callable[[], None]) -> None:
    """Run a command body; any failure renders the gapit.error/1 envelope on
    stderr and exits with the documented code (UNEXPECTED/1 for non-GapitError)."""
    try:
        action()
    except Exception as exc:
        typer.echo(render_error(exc), err=True)
        raise typer.Exit(code=exc.exit_code if isinstance(exc, GapitError) else 1) from exc
