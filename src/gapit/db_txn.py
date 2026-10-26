"""Transactional database-directory rebuilds.

A build writes several artifacts (records/sequences/index/typing/manifest).
Building in place exposes half-written databases to readers and lets two
writers corrupt each other. The fix has two layers:

1. a per-database exclusive ``flock`` serializes writers across processes;
2. the body builds in a hidden sibling staging directory, and the completed
   directory is published with a rename only after the body succeeds.

Readers therefore see the old database for the whole build, then the new one
after the swap (a short rename window remains between moving the old
directory aside and publishing the staging directory). Staging/backup
directories are hidden and ignored by database discovery.
"""

import fcntl
import os
import shutil
import tempfile
from collections.abc import Generator
from contextlib import contextmanager
from pathlib import Path

from gapit.errors import DatabaseError

# Hidden sibling prefix for lock/staging/backup entries. Discovery skips only
# this prefix, so a user database whose name happens to start with a dot is
# still visible.
TXN_PREFIX = ".gapit-txn."


def _lock_path(db_dir: Path) -> Path:
    """The hidden per-database writer lock beside the database directory."""
    return db_dir.parent / f"{TXN_PREFIX}lock.{db_dir.name}"


@contextmanager
def _database_lock(db_dir: Path) -> Generator[None, None, None]:
    """Hold the exclusive writer lock for ``db_dir`` (blocking).

    ``flock`` is released by closing the file, including on exceptions and
    process exit. The lock file itself is intentionally left in place: it is
    a stable rendezvous point, not build state.
    """
    parent = db_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    lock_path = _lock_path(db_dir)
    try:
        handle = lock_path.open("a+", encoding="utf-8")
    except OSError as exc:
        raise DatabaseError(
            f"cannot open database lock: {lock_path}",
            code="DB_LOCK_FAILED",
            context={"db": db_dir.name, "lock": str(lock_path)},
        ) from exc
    try:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        except OSError as exc:
            raise DatabaseError(
                f"cannot lock database: {lock_path}",
                code="DB_LOCK_FAILED",
                context={"db": db_dir.name, "lock": str(lock_path)},
            ) from exc
        yield
    finally:
        handle.close()


def _publish(db_dir: Path, staging: Path) -> None:
    """Swap a completed staging directory into ``db_dir``.

    The old directory is moved to a hidden backup first (same filesystem), so
    ``os.replace(staging, db_dir)`` is a rename rather than a copy. If the
    publish fails after the backup move, the old directory is restored.
    """
    backup: Path | None = None
    if db_dir.exists() or db_dir.is_symlink():
        backup = Path(
            tempfile.mkdtemp(prefix=f"{TXN_PREFIX}backup.{db_dir.name}.", dir=db_dir.parent)
        )
        backup.rmdir()
        os.replace(db_dir, backup)
    try:
        os.replace(staging, db_dir)
    except BaseException:
        if backup is not None and backup.exists():
            os.replace(backup, db_dir)
        raise
    else:
        if backup is not None:
            shutil.rmtree(backup, ignore_errors=True)


@contextmanager
def staged_database(db_dir: Path) -> Generator[Path, None, None]:
    """Build ``db_dir`` transactionally and yield the staging directory.

    The caller must write every artifact into the yielded path, never into
    ``db_dir`` directly. The canonical path is untouched while the body runs;
    on success the staging directory replaces it, on any failure (including
    ``KeyboardInterrupt``) staging is removed and the old database survives.
    Writers serialize on a hidden per-database lock, so concurrent
    ``db fetch``/``db build`` operations cannot interleave.

    The yielded directory is created with private (``0700``) permissions:
    databases are user data, and shared multi-user access should be an
    explicit ACL/chmod decision rather than the default.
    """
    parent = db_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    with _database_lock(db_dir):
        staging = Path(tempfile.mkdtemp(prefix=f"{TXN_PREFIX}build.{db_dir.name}.", dir=parent))
        try:
            yield staging
        except BaseException:
            shutil.rmtree(staging, ignore_errors=True)
            raise
        else:
            try:
                _publish(db_dir, staging)
            except BaseException:
                shutil.rmtree(staging, ignore_errors=True)
                raise
