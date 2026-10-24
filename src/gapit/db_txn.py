"""Transactional database-directory rebuilds.

A build writes several artifacts (records/sequences/index/typing/manifest).
If a ``--force`` rebuild fails after overwriting some of them, the previous
manifest can end up certifying content it never described. The standard fix
is a staging directory: move the old DB aside, build the new one at the
canonical path, and restore the old directory on any failure.
"""

import os
import shutil
import tempfile
from collections.abc import Generator
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def staged_database(db_dir: Path) -> Generator[None, None, None]:
    """Build ``db_dir`` transactionally.

    The old directory (when present) is renamed to a hidden sibling backup
    before the body runs. On success the backup is deleted; on any exception
    (including ``KeyboardInterrupt``) the partial target is removed and the
    backup is renamed back. The body must create/write into ``db_dir``.
    """
    parent = db_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    backup: Path | None = None
    if db_dir.exists():
        backup = Path(tempfile.mkdtemp(prefix=f".{db_dir.name}.backup.", dir=parent))
        backup.rmdir()
        os.replace(db_dir, backup)
    try:
        db_dir.mkdir()
        yield
    except BaseException:
        if db_dir.exists():
            shutil.rmtree(db_dir, ignore_errors=True)
        if backup is not None and backup.exists():
            os.replace(backup, db_dir)
        raise
    else:
        if backup is not None:
            shutil.rmtree(backup, ignore_errors=True)
