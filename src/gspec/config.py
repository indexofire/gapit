"""Datadir resolution and configuration defaults."""

import os
from pathlib import Path
from typing import Literal

from gspec.errors import DatabaseError

ENV_DATADIR = "GSPEC_DATADIR"
#: Pre-0.8 name, still honoured so an existing environment keeps resolving.
LEGACY_ENV_DATADIR = "GAPIT_DATADIR"
ENV_XDG_DATA_HOME = "XDG_DATA_HOME"
DEFAULT_DATADIR = Path("~/.local/share/gspec/db")
#: Pre-0.8 default. Used only when it exists and the new default does not, so
#: a 0.7 install keeps finding its databases after the rename.
LEGACY_DEFAULT_DATADIR = Path("~/.local/share/gapit/db")
DatadirSource = Literal["cli", "env", "xdg", "default"]


def _env_datadir() -> str | None:
    """The configured datadir env value: ``GSPEC_DATADIR`` then the legacy
    ``GAPIT_DATADIR``. An empty value falls through to the next name."""
    for name in (ENV_DATADIR, LEGACY_ENV_DATADIR):
        value = os.environ.get(name)
        if value:
            return value
    return None


def _default_datadir() -> Path:
    """The home default, preferring the new path but never orphaning a 0.7
    datadir: if only ``~/.local/share/gspec/db`` exists, keep using it."""
    preferred = DEFAULT_DATADIR.expanduser()
    legacy = LEGACY_DEFAULT_DATADIR.expanduser()
    if preferred.is_dir() or not legacy.is_dir():
        return preferred.resolve()
    return legacy.resolve()


def datadir_source(cli_value: Path | None) -> tuple[Path, DatadirSource]:
    """Resolve the datadir location and where the choice came from.

    Precedence: ``--datadir`` > ``$GSPEC_DATADIR`` (or the legacy
    ``$GAPIT_DATADIR``) > ``$XDG_DATA_HOME/gspec/db`` > ``~/.local/share/gspec/db``
    (falling back to a pre-0.8 ``~/.local/share/gspec/db`` when that is the
    only one present). An empty/relative
    ``$XDG_DATA_HOME`` is ignored per the XDG Base Directory specification
    (it must be an absolute path). The returned path is expanded and resolved
    but not required to exist.
    """
    if cli_value is not None:
        return Path(cli_value).expanduser().resolve(), "cli"
    env_value = _env_datadir()
    if env_value:
        return Path(env_value).expanduser().resolve(), "env"
    xdg_value = os.environ.get(ENV_XDG_DATA_HOME)
    if xdg_value:
        xdg_path = Path(xdg_value).expanduser()
        if xdg_path.is_absolute():
            return (xdg_path / "gspec" / "db").resolve(), "xdg"
    return _default_datadir(), "default"


def datadir_path(cli_value: Path | None) -> Path:
    """Resolve the datadir LOCATION without requiring it to exist.

    Precedence and XDG handling live in :func:`datadir_source`; catalog/listing
    commands use this so a fresh install can be inspected before any database
    exists. Commands that need installed databases still call
    :func:`resolve_datadir`.
    """
    return datadir_source(cli_value)[0]


def resolve_datadir(cli_value: Path | None) -> Path:
    """Resolve the datadir and require it to exist.

    Expands ``~`` and resolves to an absolute path; raises DatabaseError when
    the resolved directory does not exist. Use :func:`datadir_path` for a
    non-validating resolution (catalog listing, bootstrap diagnostics).
    """
    datadir = datadir_path(cli_value)
    if not datadir.is_dir():
        raise DatabaseError(
            f"datadir does not exist: {datadir}",
            code="DATADIR_NOT_FOUND",
            context={"datadir": str(datadir)},
        )
    return datadir


def ensure_datadir(cli_value: Path | None) -> Path:
    """Resolve the datadir and mkdir it when absent (fresh-machine bootstrap).

    Write paths (``db fetch``, ``db build``, ``db list`` catalog bootstrap,
    ``db setup``) create a missing root; every path that actually consumes an
    installed database still requires it via :func:`resolve_datadir`.
    """
    root = datadir_path(cli_value)
    try:
        # New datadirs are private by default: databases affect screening
        # results, and shared multi-user access should be an explicit ACL/chmod
        # decision rather than an accidental world-readable directory.
        root.mkdir(mode=0o700, parents=True, exist_ok=True)
    except OSError as exc:
        raise DatabaseError(
            f"cannot create datadir: {root}",
            code="DATADIR_CREATE_FAILED",
            context={"datadir": str(root)},
        ) from exc
    return root
