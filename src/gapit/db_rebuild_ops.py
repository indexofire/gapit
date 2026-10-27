"""The `db rebuild` use-case: regenerate artifacts from ``records.jsonl``.

``records.jsonl`` is the editable truth of a gene database (SPEC.md §11): the
``sequences`` FASTA, its BLAST index, and the manifest are projections of it.
Curating a database therefore means editing records — renaming a gene, fixing a
product string, merging alleles — and re-projecting. Until now that had no CLI
path of its own (``db build`` needs a FASTA input), so this use-case fills the
gap the databases guide already promised:

    Because ``records.jsonl`` is the editable truth, editing it and rebuilding
    regenerates every downstream artifact.

Scope is gene databases only. A cluster database's truth is the locus FASTA
plus ``features.json`` (derived from a source GBK/GFF), and its typing rules key
on resolved gene ids, so re-projecting from a half-edited artifact pair would
silently desynchronise them — the caller is told to rebuild from the original
GBK/GFF instead.

The rebuild is staged (:func:`gapit.db_txn.staged_database`) and reuses the
standard pipeline, so a broken records file leaves the installed database
untouched, and the manifest provenance (fetch time, source URLs, licence,
typing schema, ``source``) survives the rebuild verbatim.
"""

import shutil
from pathlib import Path
from typing import Literal

from pydantic import BaseModel

from gapit import config
from gapit.db import validate_database_name
from gapit.db_info_ops import INFO_FILENAME
from gapit.db_txn import staged_database
from gapit.dbbuild import build_database
from gapit.errors import DatabaseError
from gapit.proctools import note
from gapit.records import Manifest, count_records, read_manifest, read_records, write_records

# Sidecars a gene database may carry. typing/floors reference gene ids and the
# info sidecar holds the database's own literature/caveats; none of them carry
# the db name, so a rebuild or rename must carry them across verbatim.
SIDECARS = ("typing.json", "floors.json", INFO_FILENAME)


class RebuildReceipt(BaseModel, frozen=True):
    """One-line JSON success receipt for `db rebuild` (no biological data)."""

    db: str
    records: int
    kind: Literal["gene"]
    destination: str


def copy_sidecars(source: Path, target: Path) -> None:
    """Copy the optional typing/floors sidecars verbatim (ids, not names)."""
    for sidecar in SIDECARS:
        if (source / sidecar).is_file():
            shutil.copyfile(source / sidecar, target / sidecar)


def project_gene_records(
    source: Path,
    target: Path,
    manifest: Manifest,
    *,
    name: str | None = None,
    quiet: bool = True,
) -> int:
    """Re-project ``source/records.jsonl`` into ``target`` and build it.

    ``name`` re-points every record's ``db`` field (the ``db rename`` caller);
    None keeps the records as written (the ``db rebuild`` caller). Returns the
    record count of the rebuilt database.
    """
    db_name = name or manifest.name
    with staged_database(target) as build_dir:
        records = [
            record if name is None else record.model_copy(update={"db": db_name})
            for record in read_records(source / "records.jsonl")
        ]
        write_records(records, build_dir / "records.jsonl")
        copy_sidecars(source, build_dir)
        build_database(
            build_dir,
            name=db_name,
            dbtype=manifest.dbtype,
            source_urls=manifest.source_urls,
            fetched_at=manifest.fetched_at,
            upstream_version=manifest.upstream_version,
            content_license=manifest.license or "",
            content_note=manifest.note or "",
            typing_schema=manifest.typing_schema or "",
            source=manifest.source,
            quiet=quiet,
        )
        return count_records(build_dir / "records.jsonl")


def perform_rebuild(name: str, datadir: Path | None, *, quiet: bool = True) -> RebuildReceipt:
    """Rebuild an installed gene database in place from its own records.

    Guards: a plain name; an installed database (a manifest, so there is
    provenance to preserve); a gene kind. The rebuild is staged, so a malformed
    record or a typing reference that no longer resolves fails with the
    installed database untouched.
    """
    validate_database_name(name)
    root = config.resolve_datadir(datadir)
    db_dir = root / name
    manifest_path = db_dir / "gapit-manifest.json"
    if not manifest_path.is_file():
        raise DatabaseError(
            f"database not found: {name}",
            code="DATABASE_NOT_FOUND",
            context={"db": name, "datadir": str(root)},
        )
    manifest = read_manifest(manifest_path)
    if manifest.kind != "gene":
        raise DatabaseError(
            f"{name} is a cluster database: rebuild it from its source GBK/GFF with"
            " `gapit db build` (locus FASTA + features.json are generated artifacts,"
            " not an editable truth source)",
            code="DB_REBUILD_CLUSTER",
            context={"db": name, "kind": manifest.kind},
        )
    note(quiet, f"rebuilding {name} from records.jsonl")
    n_records = project_gene_records(db_dir, db_dir, manifest, quiet=quiet)
    return RebuildReceipt(db=name, records=n_records, kind="gene", destination=str(db_dir))
