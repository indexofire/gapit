"""Bundled database snapshots (Wave G): deterministic archives + extraction.

A snapshot is ``<name>.tar.gz`` containing EXACTLY ``records.jsonl`` and
``gapit-manifest.json`` from an installed gapit-native database directory.
card and vfdb ship inside the wheel (``src/gapit/data/snapshots/``) so the
two default databases install with zero network — snapshots eliminate the
flaky-source failure modes for exactly those two (Wave E/F3: mgc.ac.cn UA
blocks, card re-tarring, Wayback interstitials). Every other provider keeps
the upstream fetch; a future DB is bundled by dropping a ``<name>.tar.gz``
into the snapshots dir and setting ``snapshot=`` on its Provider (one line).

``records.jsonl`` (post-normalize records) is what gets snapshotted — NOT
the BLAST/minimap2 indexes: index bytes are BLAST-version-sensitive while a
local rebuild from records is deterministic and fast. The archived manifest
contributes only ``upstream_version``; fetched_at, sha256 and tool versions
are rebuilt locally by :func:`gapit.dbbuild.build_database`.

Determinism (frozen layout contract): entries added in sorted name order,
PAX format, uid/gid 0, mtime 0, and a gzip stream with mtime 0 and no
embedded filename — two builds from the same inputs are byte-identical.
"""

import gzip
import io
import os
import tarfile
from pathlib import Path
from tempfile import NamedTemporaryFile

from gapit.errors import DatabaseError
from gapit.records import Manifest

_MEMBERS = ("gapit-manifest.json", "records.jsonl")  # sorted: manifest < records


def _snapshot_error(db_dir: Path, missing: str) -> DatabaseError:
    """The uniform archive-content failure: located, machine-stable."""
    return DatabaseError(
        f"cannot snapshot {db_dir}: missing {missing}",
        code="SNAPSHOT_INVALID",
        context={"db_dir": str(db_dir), "missing": missing},
    )


def make_snapshot(db_dir: Path, dest: Path) -> None:
    """Write the deterministic snapshot archive of ``db_dir`` to ``dest``.

    Both members are required (a snapshot of a half-built db dir refuses
    loudly instead of shipping an uninstallable archive).
    """
    with dest.open("wb") as raw, gzip.GzipFile(
        # filename="" keeps the dest path OUT of the gzip header (the
        # fileobj's .name would otherwise be embedded — path-dependent bytes).
        filename="",
        mode="wb",
        fileobj=raw,
        mtime=0,
    ) as gz, tarfile.open(fileobj=gz, mode="w", format=tarfile.PAX_FORMAT) as tar:
        for name in _MEMBERS:
            source = db_dir / name
            if not source.is_file():
                raise _snapshot_error(db_dir, name)
            data = source.read_bytes()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mtime = 0
            tar.addfile(info, io.BytesIO(data))


def _member(tar: tarfile.TarFile, name: str, archive: Path) -> bytes:
    """Read one regular-file member; typed SNAPSHOT_INVALID when absent
    (no untyped KeyError escape — the Wave E card.py lesson)."""
    member = next((m for m in tar.getmembers() if m.name == name and m.isfile()), None)
    if member is None:
        raise DatabaseError(
            f"snapshot {archive} has no {name} member",
            code="SNAPSHOT_INVALID",
            context={"archive": str(archive), "missing": name},
        )
    extracted = tar.extractfile(member)
    if extracted is None:
        raise DatabaseError(
            f"{name} in {archive} is not a regular file",
            code="SNAPSHOT_INVALID",
            context={"archive": str(archive), "missing": name},
        )
    return extracted.read()


def extract_snapshot(archive: Path, db_dir: Path) -> Manifest:
    """Install ``records.jsonl`` from ``archive`` into ``db_dir`` atomically;
    return the ARCHIVED manifest (its upstream_version seeds the fresh build)."""
    with tarfile.open(archive, "r:gz") as tar:
        records = _member(tar, "records.jsonl", archive)
        archived = Manifest.model_validate_json(_member(tar, "gapit-manifest.json", archive))
    temp_path: Path | None = None
    try:
        with NamedTemporaryFile(dir=db_dir, prefix=".records.", suffix=".jsonl", delete=False) as (
            temp
        ):
            temp_path = Path(temp.name)
            temp.write(records)
        assert temp_path is not None
        os.replace(temp_path, db_dir / "records.jsonl")
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
    return archived
