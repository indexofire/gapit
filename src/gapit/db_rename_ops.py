"""The `db rename` use-case: re-home an installed database under a new name.

gapit's naming convention (genus initial + first two species letters, or the
first three genus letters) exists so database names read predictably — but the
first generation of names (`cam`, `vpa`, `cdi`, `enterotoxin`) predates it.
Renaming must not mean re-downloading or hand-editing a datadir, so this
use-case rebuilds the directory under the new name from its own truth source
and then removes the old one:

- **gene** databases carry the name inside every ``records.jsonl`` record and
  every ``gapit/v1`` header, so the directory is rebuilt through
  :func:`gapit.dbbuild.build_database` with ``db`` re-pointed (sequences,
  self-check, BLAST index, manifest — the standard pipeline, not a file copy);
- **cluster** databases carry the name *only* in the manifest (locus FASTA
  headers are bare locus ids), so their content is copied verbatim and the
  manifest is re-stamped — the committed artifacts stay byte-identical.

Both branches write through :func:`gapit.db_txn.staged_database`, so a failure
while building the new name leaves the old database fully intact; the source
directory is removed only after the target is published. Typing/floors/feature
sidecars travel along unchanged: they reference gene or locus ids, not names.
"""

import shutil
from pathlib import Path
from typing import Literal

from pydantic import BaseModel

from gapit import config
from gapit.db import make_blast_db, validate_database_name
from gapit.db_rebuild_ops import SIDECARS, project_gene_records
from gapit.db_txn import staged_database
from gapit.dbbuild import optional_tool_version_line, sha256_file, tool_version_line
from gapit.errors import DatabaseError
from gapit.proctools import note
from gapit.records import Manifest, read_manifest, write_manifest

# Sidecars a database may carry; all of them reference ids, never the db name.
_SIDECARS = (*SIDECARS, "features.json")


class RenameReceipt(BaseModel, frozen=True):
    """One-line JSON success receipt for `db rename` (no biological data)."""

    db: str
    renamed_from: str
    records: int
    kind: Literal["gene", "cluster"]
    destination: str


def _copy_sidecars(source: Path, target: Path) -> None:
    """Copy the optional typing/floors/features sidecars verbatim."""
    for sidecar in _SIDECARS:
        if (source / sidecar).is_file():
            shutil.copyfile(source / sidecar, target / sidecar)


def _rename_gene(source: Path, target: Path, manifest: Manifest, quiet: bool) -> int:
    """Re-point every gene record at the new name and rebuild the artifacts."""
    return project_gene_records(source, target, manifest, name=manifest.name, quiet=quiet)


def _rename_cluster(source: Path, target: Path, manifest: Manifest, quiet: bool) -> int:
    """Copy the locus artifacts verbatim (their headers carry no db name)."""
    with staged_database(target) as build_dir:
        shutil.copyfile(source / "sequences", build_dir / "sequences")
        _copy_sidecars(source, build_dir)
        make_blast_db(build_dir / "sequences", manifest.name, dbtype=manifest.dbtype)
        note(quiet, f"BLAST index built ({manifest.dbtype})")
        write_manifest(
            Manifest(
                name=manifest.name,
                source_urls=manifest.source_urls,
                fetched_at=manifest.fetched_at,
                sha256=sha256_file(build_dir / "sequences"),
                n_records=manifest.n_records,
                dbtype=manifest.dbtype,
                kind="cluster",
                header_format="plain",
                upstream_version=manifest.upstream_version,
                typing_schema=manifest.typing_schema,
                license=manifest.license,
                note=manifest.note,
                source=manifest.source,
                makeblastdb_version=tool_version_line(["makeblastdb", "-version"]),
                minimap2_version=optional_tool_version_line(["minimap2", "--version"]),
            ),
            build_dir / "gapit-manifest.json",
        )
        return manifest.n_records


def perform_rename(
    old: str,
    new: str,
    datadir: Path | None,
    *,
    force: bool = False,
    quiet: bool = True,
) -> RenameReceipt:
    """Rename an installed database ``old`` -> ``new`` under the datadir.

    Guards, in order: both names must be plain names (no path separators or
    control characters); they must differ; the source must be an installed
    gapit database (its manifest supplies the provenance the rebuild keeps —
    an abricate-style directory has no name to rewrite); the target must not
    already exist unless ``force``.
    """
    validate_database_name(old)
    validate_database_name(new)
    if old == new:
        raise DatabaseError(
            f"database is already named {new!r}",
            code="DB_ALREADY_EXISTS",
            context={"db": new},
        )
    root = config.resolve_datadir(datadir)
    source = root / old
    target = root / new
    manifest_path = source / "gapit-manifest.json"
    if not manifest_path.is_file():
        raise DatabaseError(
            f"database not found: {old}",
            code="DATABASE_NOT_FOUND",
            context={"db": old, "datadir": str(root)},
        )
    if target.exists() and not force:
        raise DatabaseError(
            f"won't overwrite existing database {new} (use --force)",
            code="DB_ALREADY_EXISTS",
            context={"db": new},
        )
    manifest = read_manifest(manifest_path).model_copy(update={"name": new})
    note(quiet, f"renaming {old} -> {new} ({manifest.kind} database)")
    rename = _rename_cluster if manifest.kind == "cluster" else _rename_gene
    n_records = rename(source, target, manifest, quiet)
    try:
        shutil.rmtree(source)
    except OSError as exc:
        raise DatabaseError(
            f"renamed to {new} but could not remove the old directory {source}: {exc}",
            code="FILE_IO_ERROR",
            context={"db": new, "old": str(source)},
        ) from exc
    return RenameReceipt(
        db=new,
        renamed_from=old,
        records=n_records,
        kind=manifest.kind,
        destination=str(target),
    )
