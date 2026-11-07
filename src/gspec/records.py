"""The per-database truth source: Record JSONL store and gapit.manifest/1.

``records.jsonl`` is the editable truth a database directory is built from
(the ``sequences`` FASTA and its BLAST index are generated projections, built
later in Wave A3); ``gspec-manifest.json`` is the provenance sidecar describing
how that truth was obtained. Both are FILE contracts — they never appear on
stdout and are not registered in ``gspec schema``.
"""

import re
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator

from gspec.db_txn import TXN_PREFIX
from gspec.errors import InputError
from gspec.formats.json import ToolDocument
from gspec.layout import manifest_path
from gspec.schema_ids import Canonical

_SHA256_SHAPE = re.compile(r"[0-9a-fA-F]{64}")


class Record(BaseModel, frozen=True):
    """One known gene — the unit of truth a database projection is built from.

    ``function`` carries functional categories (AMR classes, virulence,
    O-antigen, replicon, biocide) — gspec answers presence/absence, so the
    slot's old ``resistance`` name was a misnomer (Wave F1).
    """

    db: str
    gene: str
    sequence: str
    accession: str = ""
    function: tuple[str, ...] = ()
    product: str = "n/a"
    source_id: str = ""


def write_records(records: Iterable[Record], path: Path) -> None:
    """Write records as deterministic JSONL: given order, one JSON object per
    line, LF endings, UTF-8."""
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(record.model_dump_json())
            handle.write("\n")


def read_records(path: Path) -> Iterator[Record]:
    """Stream records back line by line (nothing is materialized).

    Blank lines are skipped; a non-parsing line raises InputError
    RECORDS_MALFORMED with the physical line number in context. The path is
    checked eagerly; iteration itself is lazy.
    """
    if not path.is_file():
        raise InputError(
            f"records file not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        )
    return _iter_records(path)


def _iter_records(path: Path) -> Iterator[Record]:
    try:
        handle = path.open(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise InputError(
            f"records file not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        ) from exc
    with handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                yield Record.model_validate_json(line)
            except ValidationError as exc:
                raise InputError(
                    f"malformed record at {path}:{line_number}: {_reason(exc)}",
                    code="RECORDS_MALFORMED",
                    context={"file": str(path), "line": str(line_number)},
                ) from exc


def count_records(path: Path) -> int:
    """Count record lines by streaming (blank lines excluded, no list built)."""
    return sum(1 for _ in read_records(path))


class Manifest(BaseModel, frozen=True):
    """gapit.manifest/1 — provenance sidecar for one database directory."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Annotated[Literal["gspec.manifest/1"], Canonical] = Field(
        default="gspec.manifest/1", alias="schema"
    )
    name: str
    source_urls: tuple[str, ...]
    fetched_at: str
    sha256: str
    n_records: int = Field(ge=0)
    dbtype: Literal["nucl", "prot"]
    # gene: per-gene records pipeline (default; legacy dirs are gene dbs).
    # cluster: GBK/GFF-built gene-cluster db (features.json + locus FASTA).
    kind: Literal["gene", "cluster"] = "gene"
    # gspec/v1: tagged per-gene headers; plain: bare locus-id headers (cluster).
    # gapit/v1 is the pre-0.9 spelling and is still accepted on read.
    header_format: Literal["gspec/v1", "plain"] = "gspec/v1"

    @field_validator("header_format", mode="before")
    @classmethod
    def _legacy_header_format(cls, value: object) -> object:
        return "gspec/v1" if value == "gapit/v1" else value

    upstream_version: str = ""
    # Typing provenance: the gspec.typing schema tag of the installed
    # typing.json (None when the database is untyped; omitted at
    # serialization so pre-typing manifests stay byte-identical).
    typing_schema: str | None = None
    # Database-content provenance (kaptive-style cluster providers); None is
    # omitted at serialization so pre-license manifests stay byte-identical.
    license: str | None = None
    note: str | None = None
    # "bundled" marks a database materialized from the wheel's own content
    # (gspec.bundled — public-domain panels); None is omitted at
    # serialization so every other manifest stays byte-identical.
    source: Literal["bundled"] | None = None
    tool: ToolDocument = ToolDocument()
    makeblastdb_version: str = ""
    minimap2_version: str = ""

    @field_validator("sha256")
    @classmethod
    def _sha256_hex(cls, value: str) -> str:
        """Exactly 64 hex chars; case-insensitive on input, stored lowercase."""
        if not _SHA256_SHAPE.fullmatch(value):
            raise ValueError("sha256 must be exactly 64 hex characters")
        return value.lower()


def _reason(exc: ValidationError) -> str:
    """The validation failure rendered as one line (str(exc) is multiline)."""
    return " ".join(str(exc).split())


def write_manifest(manifest: Manifest, path: Path) -> None:
    """Write the manifest as indented JSON (by alias, ``schema`` first) with a
    trailing newline, LF endings, UTF-8. None-valued optional fields
    (license/note) are omitted — every pre-existing manifest byte-identical."""
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(manifest.model_dump_json(indent=2, by_alias=True, exclude_none=True))
        handle.write("\n")


def read_manifest(path: Path) -> Manifest:
    """Parse and validate a gapit.manifest/1 file."""
    if not path.is_file():
        raise InputError(
            f"manifest file not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        )
    try:
        return Manifest.model_validate_json(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError) as exc:
        raise InputError(
            f"manifest file not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        ) from exc
    except ValidationError as exc:
        raise InputError(
            f"malformed manifest {path}: {_reason(exc)}",
            code="MANIFEST_MALFORMED",
            context={"file": str(path)},
        ) from exc


def installed_db_dirs(root: Path) -> list[Path]:
    """Datadir subdirectories holding a manifest, sorted by name."""
    return sorted(
        (
            child
            for child in root.iterdir()
            if not child.name.startswith(TXN_PREFIX)
            and child.is_dir()
            and manifest_path(child).is_file()
        ),
        key=lambda child: child.name,
    )
