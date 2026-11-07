"""Database-directory layout names shared by the build, read and list paths.

Kept dependency-free (pathlib only) so both :mod:`gspec.db` and
:mod:`gspec.records` can import it without a cycle.
"""

from __future__ import annotations

from pathlib import Path

#: Manifest filename written by 0.9 and later.
MANIFEST_NAME = "gspec-manifest.json"

#: Pre-0.9 filename, still read so an existing datadir keeps working.
LEGACY_MANIFEST_NAME = "gapit-manifest.json"


def manifest_path(db_dir: Path) -> Path:
    """The manifest file for a database directory.

    Prefers the current name, falls back to the pre-0.9 one. A directory with
    neither returns the current name so callers can test ``is_file()``.
    """
    current = db_dir / MANIFEST_NAME
    if current.is_file():
        return current
    legacy = db_dir / LEGACY_MANIFEST_NAME
    return legacy if legacy.is_file() else current


__all__ = ["LEGACY_MANIFEST_NAME", "MANIFEST_NAME", "manifest_path"]
