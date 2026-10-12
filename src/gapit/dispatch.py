"""Shared CLI plumbing: error-envelope dispatch, the ``--datadir`` and
``--quiet`` option factories, the Enum-backed choice type, and the
``--output`` delivery target (stdout chunks vs a lazily opened file).

Imported by every command module (cli.py, cmd_*.py); imports nothing from
gapit except :mod:`gapit.errors`, so it can never participate in a cycle.
"""

import enum
from collections.abc import Callable, Generator
from contextlib import contextmanager
from pathlib import Path
from typing import IO, Any, TextIO, TypeVar

import click

from gapit.errors import GapitError, render_error

# The option factories return click's pass-through option decorator; the
# Any-bound TypeVar mirrors click.decorators' own FC typing (decorated
# functions keep their concrete signatures).
Decorated = TypeVar("Decorated", bound=Callable[..., Any])

EnumT = TypeVar("EnumT", bound=enum.Enum)


class EnumChoice(click.Choice[EnumT]):
    """Choice over an Enum's VALUES that still returns the member.

    click's plain Choice normalizes Enum choices by member NAME (``map_ont``),
    which would reject hyphenated CLI tokens like ``--read-type map-ont``;
    normalizing on ``value`` keeps the token set frozen and hands back the
    member the function signatures are annotated with.
    """

    def __init__(self, choices: type[EnumT]) -> None:
        super().__init__(list(choices))

    def normalize_choice(self, choice: object, ctx: click.Context | None) -> str:
        return str(getattr(choice, "value", choice))


def datadir_option() -> Callable[[Decorated], Decorated]:
    """The shared ``--datadir`` option (one help text across every command)."""
    return click.option(
        "--datadir",
        "-D",
        "datadir",
        type=click.Path(path_type=Path),
        help="Database directory (default: $GAPIT_DATADIR, then ~/.local/share/gapit/db).",
    )


def quiet_option() -> Callable[[Decorated], Decorated]:
    """The shared ``--quiet`` flag (byte-identical across commands)."""
    return click.option(
        "--quiet",
        "-q",
        is_flag=True,
        help="Silence stderr diagnostics.",
    )


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
        yield lambda chunk: click.echo(chunk, nl=False)
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
        click.echo(render_error(exc), err=True)
        raise SystemExit(exc.exit_code if isinstance(exc, GapitError) else 1) from exc
