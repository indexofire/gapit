"""Datadir resolution and configuration defaults."""

import os
from pathlib import Path

from gapit.errors import DatabaseError

ENV_DATADIR = "GAPIT_DATADIR"
DEFAULT_DATADIR = Path("~/.local/share/gapit/db")


def resolve_datadir(cli_value: Path | None) -> Path:
    """Resolve the datadir: CLI ``--datadir`` > ``$GAPIT_DATADIR`` > ``~/.local/share/gapit/db``.

    Expands ``~`` and resolves to an absolute path; raises DatabaseError when the
    resolved directory does not exist.
    """
    raw = cli_value if cli_value is not None else os.environ.get(ENV_DATADIR, DEFAULT_DATADIR)
    datadir = Path(raw).expanduser().resolve()
    if not datadir.is_dir():
        raise DatabaseError(
            f"datadir does not exist: {datadir}",
            code="DATADIR_NOT_FOUND",
            context={"datadir": str(datadir)},
        )
    return datadir


def ensure_datadir(cli_value: Path | None) -> Path:
    """Resolve the datadir and mkdir it when absent (fresh-machine bootstrap).

    Write paths (``db fetch``, ``db build``) bootstrap a missing root; every
    read path still demands it via ``resolve_datadir``. The resolved path is
    recovered from ``resolve_datadir``'s DATADIR_NOT_FOUND context — this
    module stays the single owner of resolution.
    """
    try:
        root = resolve_datadir(cli_value)
    except DatabaseError as exc:
        if exc.code != "DATADIR_NOT_FOUND":
            raise
        root = Path(exc.context["datadir"])
    try:
        root.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise DatabaseError(
            f"cannot create datadir: {root}",
            code="DATADIR_CREATE_FAILED",
            context={"datadir": str(root)},
        ) from exc
    return root
