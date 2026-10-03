"""Shared stderr-bound rich progress bar for long batch runs (gmlst design).

Modeled on gmlst's ``make_progress``: spinner, description, bar, M-of-N,
percent, elapsed time — on a console bound to stderr so bar rendering can
never pollute stdout data streams (gapit's stdout-purity law). One
deliberate deviation from the gmlst constructor: rich's Live normally
REDIRECTS the process's real stdout/stderr through the bar's console on a
terminal, which would reroute streamed report bytes onto stderr here — so
``redirect_stdout=False`` / ``redirect_stderr=False`` are pinned.

The bar is ACTIVE only when ``--quiet`` is off and stderr is a TTY. When
inactive, :func:`screen_progress` yields a no-op shim (task id ``None``)
so call sites keep one code path and the legacy per-file stderr notes
print byte-identically; when active the bar REPLACES those notes (its
M/N + description carry them).
"""

from collections.abc import Generator, Sequence
from contextlib import contextmanager
from typing import Final, Protocol

from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TaskID,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)

# Module-level so tests can substitute a recording (record=True) console.
STATUS_CONSOLE: Final = Console(stderr=True)


class BatchProgress(Protocol):
    """The single-task bar slice the batch call sites use."""

    def describe(self, text: str) -> None:
        """Set the bar's description (e.g. the file/sample/db being worked)."""
        ...

    def advance(self) -> None:
        """Mark one batch item completed."""
        ...


class IdleBatch:
    """No-op shim used when the bar is disabled (quiet or non-TTY stderr).

    Absorbs describe/advance so call sites stay branch-free; the legacy
    per-file stderr notes keep printing at the call sites exactly as
    before.
    """

    def describe(self, text: str) -> None:
        return None

    def advance(self) -> None:
        return None


class _LiveBatch:
    """Adapter pinning a rich Progress to its one batch task."""

    def __init__(self, progress: Progress, task_id: TaskID) -> None:
        self._progress = progress
        self._task_id = task_id

    def describe(self, text: str) -> None:
        self._progress.update(self._task_id, description=text)

    def advance(self) -> None:
        self._progress.advance(self._task_id)


_IDLE: Final = IdleBatch()


def make_progress() -> Progress:
    """The gmlst column set (spinner, description, bar, M/N, %, elapsed) on
    the stderr console; stdout/stderr redirection stays OFF (see module
    docstring) so the process's stdout data stream is never captured."""
    return Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        console=STATUS_CONSOLE,
        redirect_stdout=False,
        redirect_stderr=False,
    )


@contextmanager
def screen_progress(
    items: Sequence[object], *, quiet: bool, description: str = "Screening"
) -> Generator[tuple[BatchProgress, TaskID | None], None, None]:
    """Drive a batch over ``items`` behind the bar; yields ``(bar,
    task_id)``.

    ``task_id`` is None when the bar is INACTIVE — under ``--quiet`` or
    when stderr is not a TTY (piped/CI output keeps logs clean) — and then
    ``bar`` is the no-op shim: call sites describe/advance freely and keep
    printing the legacy per-file notes. Active, the task starts at
    ``description`` with total ``len(items)`` and call sites REPLACE their
    per-file notes with bar updates.
    """
    progress = None if quiet else make_progress()
    if progress is not None and progress.console.is_terminal:
        with progress:
            task_id = progress.add_task(description, total=len(items))
            yield _LiveBatch(progress, task_id), task_id
        return
    yield _IDLE, None
