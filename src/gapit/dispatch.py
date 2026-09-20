"""Shared CLI plumbing: error-envelope dispatch and the ``--datadir`` option.

Imported by every command module (cli.py, cmd_*.py); imports nothing from
gapit except :mod:`gapit.errors`, so it can never participate in a cycle.
"""

from collections.abc import Callable
from pathlib import Path
from typing import Annotated

import typer

from gapit.errors import GapitError, render_error

Datadir = Annotated[
    Path | None,
    typer.Option(
        "--datadir",
        help="Database directory (default: $GAPIT_DATADIR, then ~/.local/share/gapit/db).",
    ),
]


def dispatch(action: Callable[[], None]) -> None:
    """Run a command body; any failure renders the gapit.error/1 envelope on
    stderr and exits with the documented code (UNEXPECTED/1 for non-GapitError)."""
    try:
        action()
    except Exception as exc:
        typer.echo(render_error(exc), err=True)
        raise typer.Exit(code=exc.exit_code if isinstance(exc, GapitError) else 1) from exc
