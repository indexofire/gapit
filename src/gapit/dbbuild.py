"""Wave A3 build pipeline: records.jsonl -> a fully built gapit-native database.

``build_database`` turns a directory's truth source (``records.jsonl``) into
the three artifacts consumers rely on — the ``sequences`` FASTA projection,
its BLAST index, and the ``gapit-manifest.json`` provenance sidecar written
LAST (it certifies the artifacts). Every step is deterministic, atomic where
it matters, and self-verifying: the generated headers must decode back
through :mod:`gapit.dbcodec` before anything is indexed. No minimap2 index
is persisted: reads mode indexes the ``sequences`` FASTA in memory only
(SPEC.md §10).
"""

import hashlib
import os
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Literal

from gapit.db import DbHeader, make_blast_db
from gapit.dbcodec import decode_seqid, encode_seqid
from gapit.errors import DatabaseError, DependencyError
from gapit.fasta import iter_fasta
from gapit.records import Manifest, Record, count_records, read_records, write_manifest

_WRAP_COLUMNS = 60
_CHUNK_BYTES = 1 << 20


def _run(argv: list[str]) -> subprocess.CompletedProcess[str]:
    """Run an external tool with an argv list (never a shell).

    Mirrors gapit.db._run: db.py is frozen in this wave except for the
    make_blast_db dbtype parameter, so the runner lives here too.
    """
    try:
        return subprocess.run(argv, check=False, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise DependencyError(
            f"required binary not found on PATH: {argv[0]}",
            code="MISSING_DEPENDENCY",
            context={"binary": argv[0]},
        ) from exc


def generate_sequences(records_path: Path, sequences_path: Path) -> None:
    """Stream ``records.jsonl`` into the ``sequences`` FASTA projection.

    One header line per record (``encode_seqid(...)`` + space + product),
    sequence wrapped at 60 columns, LF endings, UTF-8. Bytes land in a temp
    file in the target directory and are ``os.replace``d only after the whole
    file was written, so a failure never leaves a partial ``sequences``.

    Raises DatabaseError ``BUILD_INVALID`` (exit 4) for a record whose product
    contains a line break or whose sequence is empty, and for a records file
    with zero records (makeblastdb cannot index nothing); the per-record
    errors carry the record's gene in context. Duplicate ``(db, gene)`` pairs
    are allowed and preserved — upstream databases contain them and
    records.jsonl order is the truth.
    """
    temp_path: Path | None = None
    n_written = 0
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
            for record in read_records(records_path):
                _check_record(record)
                seqid = encode_seqid(record.db, record.gene, record.accession, record.function)
                temp.write(f">{seqid} {record.product}\n")
                for offset in range(0, len(record.sequence), _WRAP_COLUMNS):
                    temp.write(f"{record.sequence[offset : offset + _WRAP_COLUMNS]}\n")
                n_written += 1
        if n_written == 0:
            raise DatabaseError(
                f"cannot build {sequences_path}: {records_path} contains no records",
                code="BUILD_INVALID",
                context={"file": str(records_path)},
            )
        assert temp_path is not None
        os.replace(temp_path, sequences_path)
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)


def _check_record(record: Record) -> None:
    """Generation-time validation: the two record shapes that would produce a
    broken FASTA (headers must stay one line; empty sequences cannot be
    indexed)."""
    if "\n" in record.product or "\r" in record.product:
        raise DatabaseError(
            f"record {record.gene!r}: product contains a line break",
            code="BUILD_INVALID",
            context={"gene": record.gene},
        )
    if not record.sequence:
        raise DatabaseError(
            f"record {record.gene!r}: empty sequence",
            code="BUILD_INVALID",
            context={"gene": record.gene},
        )


def _self_check_failed(gene: str, reason: str) -> DatabaseError:
    """The uniform self-check failure: gene located, machine-stable reason."""
    return DatabaseError(
        f"self-check failed at gene {gene!r} ({reason})",
        code="BUILD_SELF_CHECK_FAILED",
        context={"gene": gene, "reason": reason},
    )


def verify_sequences(sequences_path: Path, records_path: Path, name: str) -> None:
    """Self-check: every generated FASTA record decodes back to its record.

    Streams ``sequences`` and ``records.jsonl`` in parallel; the header (via
    :func:`gapit.dbcodec.decode_seqid` with ``default_db=name``) and the
    sequence must match record-for-record, in order. The FASTA description
    (product) is deliberately not compared: the reader strips leading and
    trailing whitespace by design, so it is a display field, not a truth
    field. Any decode failure, field mismatch, or count mismatch raises
    ``BUILD_SELF_CHECK_FAILED`` (a structurally invalid FASTA surfaces as its
    native InputError instead — either way the build stops before indexing).
    """
    fasta_stream = iter_fasta(sequences_path)
    for record in read_records(records_path):
        fasta = next(fasta_stream, None)
        if fasta is None:
            raise _self_check_failed(record.gene, "missing_fasta_record")
        try:
            header = decode_seqid(fasta.id, default_db=name)
        except DatabaseError as exc:
            raise _self_check_failed(record.gene, f"decode:{exc.code}") from exc
        expected = DbHeader(
            database=record.db or name,
            gene=record.gene,
            accession=record.accession,
            function=";".join(record.function),
        )
        if header != expected:
            raise _self_check_failed(record.gene, "header_mismatch")
        if fasta.sequence != record.sequence:
            raise _self_check_failed(record.gene, "sequence_mismatch")
    if next(fasta_stream, None) is not None:
        raise _self_check_failed("", "extra_fasta_record")


def _sha256(path: Path) -> str:
    """Streaming SHA256 of a file's bytes (1 MiB chunks)."""
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(_CHUNK_BYTES):
            hasher.update(chunk)
    return hasher.hexdigest()


def _version_line(argv: list[str]) -> str:
    """First line of a version command's stdout ('' when it printed nothing)."""
    stdout = _run(argv).stdout
    return stdout.splitlines()[0].strip() if stdout else ""


def _note(quiet: bool, message: str) -> None:
    """Per-step progress on stderr when quiet is disabled (stdout stays pure)."""
    if not quiet:
        print(f"gapit: {message}", file=sys.stderr)


def build_database(
    db_dir: Path,
    *,
    name: str,
    dbtype: Literal["nucl", "prot"],
    source_urls: Sequence[str],
    fetched_at: str,
    upstream_version: str = "",
    quiet: bool = True,
    debug: bool = False,
) -> Manifest:
    """Build every gapit-native artifact in ``db_dir`` from its records.jsonl.

    Pipeline — each step runs only after the previous one succeeded, and the
    manifest is written LAST (it certifies the artifacts):

    1. ``generate_sequences`` -> ``db_dir/sequences``
    2. self-check: every header decodes back (``BUILD_SELF_CHECK_FAILED``)
    3. streaming SHA256 of ``sequences``
    4. ``makeblastdb`` with the EXPLICIT ``dbtype`` — the manifest declares
       the type, so the mol_type heuristic is skipped
    5. count records and capture ``blastn -version`` / ``minimap2 --version``
    6. write ``db_dir/gapit-manifest.json`` and return the Manifest

    No ``.mmi`` is built for either dbtype: reads mode indexes the FASTA in
    memory with the invocation preset's own parameters (minimap2 is
    nucleotide-only, which is why prot databases never participated anyway).
    """
    records_path = db_dir / "records.jsonl"
    sequences_path = db_dir / "sequences"
    generate_sequences(records_path, sequences_path)
    _note(quiet, f"generated {sequences_path}")
    verify_sequences(sequences_path, records_path, name)
    _note(quiet, f"self-check passed for {name}")
    sha256 = _sha256(sequences_path)
    make_blast_db(sequences_path, name, dbtype=dbtype, debug=debug)
    _note(quiet, f"BLAST index built ({dbtype})")
    manifest = Manifest(
        name=name,
        source_urls=tuple(source_urls),
        fetched_at=fetched_at,
        sha256=sha256,
        n_records=count_records(records_path),
        dbtype=dbtype,
        upstream_version=upstream_version,
        makeblastdb_version=_version_line(["blastn", "-version"]),
        # Kept although no .mmi is built: environment provenance for the
        # machine that produced the artifacts (spec'd in gapit.manifest/1).
        minimap2_version=_version_line(["minimap2", "--version"]),
    )
    write_manifest(manifest, db_dir / "gapit-manifest.json")
    return manifest
