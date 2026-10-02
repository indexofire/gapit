"""Positional-input confusion guard for the ``gapit screen`` command.

Lives beside :mod:`gapit.cmd_screen` (same layer as
:mod:`gapit.cmd_screen_reads_args`): catches the word-style slip where a
database NAME lands in the positionals, so the failure names ``--db``/``-d``
instead of the generic input error plus a starved downstream pipe.
"""

from pathlib import Path

from gapit.bundled import bundled_names
from gapit.config import resolve_datadir
from gapit.db import discover_databases
from gapit.errors import DatabaseError, UsageError
from gapit.providers import REGISTRY


def reject_db_name_positionals(files: list[Path] | None, datadir: Path | None) -> None:
    """Fail with a targeted usage error when a positional FILE is missing on
    disk but names a known database (provider catalog, bundled wheel content,
    or datadir install) — e.g. click eating ``-db ecoli_dec`` as ``-d b`` plus
    a positional. Positionals that exist, or missing files no database
    answers to, pass through untouched (ordinary INPUT_NOT_FOUND)."""
    missing = [path for path in files or [] if not path.is_file()]
    if not missing:
        return
    names = set(REGISTRY) | set(bundled_names())
    try:
        names |= {database.name for database in discover_databases(resolve_datadir(datadir))}
    except DatabaseError as exc:
        if exc.code != "DATADIR_NOT_FOUND":
            raise
        # a missing datadir still leaves catalog + bundled names; the use-case
        # raises the typed datadir error in its own order
    for path in missing:
        if path.name in names:
            raise UsageError(
                f"input file not found: '{path.name}' — it is a database NAME; screen takes"
                f" databases via --db/-d/-db (e.g. --db {path.name}), positional arguments"
                " are genome files",
                code="USAGE_ERROR",
                context={"file": str(path), "database": path.name},
            )
