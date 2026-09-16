"""Datadir resolution and configuration defaults."""

import os
from pathlib import Path

from gapit.errors import DatabaseError

ENV_DATADIR = "GAITA_DATADIR"
DEFAULT_DATADIR = Path("~/.local/share/gapit/db")


def resolve_datadir(cli_value: Path | None) -> Path:
    """Resolve the datadir: CLI ``--datadir`` > ``$GAITA_DATADIR`` > ``~/.local/share/gapit/db``.

    Expands ``~`` and resolves to an absolute path; raises DatabaseError when the
    resolved directory does not exist.
    """
    raw = cli_value if cli_value is not None else os.environ.get(ENV_DATADIR, DEFAULT_DATADIR)
    datadir = Path(raw).expanduser().resolve()
    if not datadir.is_dir():
        raise DatabaseError(f"datadir does not exist: {datadir}", code="DATADIR_NOT_FOUND")
    return datadir
