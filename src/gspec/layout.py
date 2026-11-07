"""Database-directory layout names shared by the build, read and list paths.

Kept dependency-free (pathlib only) so both :mod:`gspec.db` and
:mod:`gspec.records` can import it without a cycle.
"""

from __future__ import annotations

from pathlib import Path

#: Manifest filename written by 0.9 and later.
MANIFEST_NAME = "gspec-manifest.json"


def manifest_path(db_dir: Path) -> Path:
    """The manifest file of a database directory."""
    return db_dir / MANIFEST_NAME


__all__ = ["MANIFEST_NAME", "manifest_path"]
