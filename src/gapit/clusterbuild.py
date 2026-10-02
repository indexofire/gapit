"""The cluster-db build pipeline with self-check (the gapit.dbbuild sibling).

Turns parsed loci into a cluster database directory: locus FASTA
``sequences`` (one bare-id record per locus), the ``gapit.features/1``
feature table (``features.json``), an optional validated ``gapit.typing``
spec (v1 or v2; cluster dbs evaluate one scheme) copied in as
``typing.json``, the makeblastdb nucl index (the gene-view
screening surface for stage 2), and the manifest — written LAST, certifying
every artifact, with ``kind: cluster`` and plain locus-id headers.

Input kind is suffix-based: GBK family (.gbk/.gbff/.gb[.gz/.bz2]) parses
through gapit.gbfeatures, GFF family (.gff/.gff3[.gz/.bz2]) through
gapit.gffparse; every other suffix is a gene-pipeline input (db_build_ops
routes those to gapit.dbbuild). ``--typing`` specs are validated BEFORE any
artifact is written, so a malformed spec leaves no half-built database.
"""

import os
import shutil
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path
from tempfile import NamedTemporaryFile

from gapit.db import make_blast_db
from gapit.dbbuild import sha256_file, tool_version_line
from gapit.errors import DatabaseError
from gapit.fasta import iter_fasta
from gapit.gbfeatures import (
    FeaturesDocument,
    LocusSequence,
    parse_genbank_features,
    write_features_document,
)
from gapit.gffparse import parse_gff3_features
from gapit.proctools import note
from gapit.records import Manifest, write_manifest
from gapit.typing_models import (
    read_typing_document,
    single_scheme,
    typing_schema_of,
    validate_references,
)

_WRAP_COLUMNS = 60


def _with_compression(suffixes: tuple[str, ...]) -> tuple[str, ...]:
    """A suffix tuple plus its .gz/.bz2 wrapped variants."""
    return suffixes + tuple(f"{s}{ext}" for ext in (".gz", ".bz2") for s in suffixes)


_GBK_SUFFIXES: tuple[str, ...] = (".gbk", ".gbff", ".gb")
_GFF_SUFFIXES: tuple[str, ...] = (".gff", ".gff3")
_COMPRESSED_SUFFIXES = _with_compression(_GBK_SUFFIXES + _GFF_SUFFIXES)


def is_cluster_input(path: Path) -> bool:
    """Suffix sniff: GBK/GFF-family inputs (plain or compressed) build
    cluster databases; everything else runs the gene pipeline."""
    return path.name.lower().endswith(_COMPRESSED_SUFFIXES)


def parse_cluster_input(path: Path) -> list[LocusSequence]:
    """Dispatch the cluster parser by suffix (GFF family vs GBK family)."""
    if path.name.lower().endswith(_with_compression(_GFF_SUFFIXES)):
        return parse_gff3_features(path)
    return parse_genbank_features(path)


def _check_failed(locus: str, reason: str) -> DatabaseError:
    """The cluster-build self-check failure: locus located, stable reason."""
    return DatabaseError(
        f"self-check failed at locus {locus!r} ({reason})",
        code="BUILD_SELF_CHECK_FAILED",
        context={"locus": locus, "reason": reason},
    )


def _write_locus_fasta(loci: Sequence[LocusSequence], sequences_path: Path) -> None:
    """Stream parsed loci into the ``sequences`` FASTA: one bare-id header
    per locus, sequence wrapped at 60 columns (the generate_sequences
    atomicity pattern — temp file + os.replace, no partial ``sequences``)."""
    temp_path: Path | None = None
    try:
        with NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            dir=sequences_path.parent,
            prefix=f".{sequences_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as temp:
            temp_path = Path(temp.name)
            for locus in loci:
                temp.write(f">{locus.features.id}\n")
                for offset in range(0, len(locus.sequence), _WRAP_COLUMNS):
                    temp.write(f"{locus.sequence[offset : offset + _WRAP_COLUMNS]}\n")
        assert temp_path is not None
        os.replace(temp_path, sequences_path)
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)


def _verify_cluster(loci: Sequence[LocusSequence], sequences_path: Path) -> None:
    """Self-check: every locus round-trips through the written FASTA — same
    count, same order, same id, same sequence bytes (gene bounds were proven
    at parse time; this certifies the artifact)."""
    fasta_stream = iter_fasta(sequences_path)
    for locus in loci:
        fasta = next(fasta_stream, None)
        if fasta is None:
            raise _check_failed(locus.features.id, "missing_fasta_record")
        if fasta.id != locus.features.id:
            raise _check_failed(locus.features.id, "id_mismatch")
        if fasta.sequence != locus.sequence:
            raise _check_failed(locus.features.id, "sequence_mismatch")
    if next(fasta_stream, None) is not None:
        raise _check_failed("", "extra_fasta_record")


def build_cluster_database(
    db_dir: Path,
    *,
    name: str,
    loci: Sequence[LocusSequence],
    typing: Path | None,
    source_urls: Sequence[str],
    fetched_at: str,
    typing_schema: str = "",
    quiet: bool = True,
    debug: bool = False,
    content_license: str = "",
    content_note: str = "",
) -> Manifest:
    """Build every cluster-db artifact in ``db_dir`` (manifest written LAST):
    locus FASTA -> features.json -> self-check -> optional typing copy ->
    makeblastdb (nucl; cluster loci are always nucleotide) -> manifest.
    ``typing_schema`` records the installed spec's gapit.typing version in
    the manifest. ``content_license``/``content_note`` carry database-content
    provenance (kaptive-style providers); both default to empty and
    serialize only when set."""
    sequences_path = db_dir / "sequences"
    _write_locus_fasta(loci, sequences_path)
    note(quiet, f"generated {sequences_path}")
    write_features_document(tuple(locus.features for locus in loci), db_dir / "features.json")
    _verify_cluster(loci, sequences_path)
    note(quiet, f"self-check passed for {name}")
    if typing is not None:
        shutil.copyfile(typing, db_dir / "typing.json")
    sha256 = sha256_file(sequences_path)
    make_blast_db(sequences_path, name, dbtype="nucl", debug=debug)
    note(quiet, "BLAST index built (nucl)")
    manifest = Manifest(
        name=name,
        source_urls=tuple(source_urls),
        fetched_at=fetched_at,
        sha256=sha256,
        n_records=len(loci),
        dbtype="nucl",
        kind="cluster",
        header_format="plain",
        typing_schema=typing_schema or None,
        makeblastdb_version=tool_version_line(["blastn", "-version"]),
        minimap2_version=tool_version_line(["minimap2", "--version"]),
        license=content_license or None,
        note=content_note or None,
    )
    write_manifest(manifest, db_dir / "gapit-manifest.json")
    return manifest


def perform_cluster_build(
    name: str, fasta: Path, typing: Path | None, db_dir: Path, *, quiet: bool = True
) -> Manifest:
    """The cluster branch of the build use-case: parse the GBK/GFF input,
    validate the typing spec (schema, referenced ids, AND exactly one
    scheme — gapit.cluster/1 has one phenotype slot — so a bad spec fails
    before any artifact is written), then build (callers shape the
    receipt)."""
    loci = parse_cluster_input(fasta)
    typing_schema = ""
    if typing is not None:
        document = read_typing_document(typing)
        single_scheme(document)
        validate_references(document, FeaturesDocument(loci=tuple(item.features for item in loci)))
        typing_schema = typing_schema_of(typing)
    db_dir.mkdir(parents=True, exist_ok=True)
    return build_cluster_database(
        db_dir,
        name=name,
        loci=loci,
        typing=typing,
        source_urls=("local",),
        fetched_at=datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        typing_schema=typing_schema,
        quiet=quiet,
    )
