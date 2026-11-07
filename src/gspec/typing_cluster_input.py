"""`gspec typing` input for cluster screen results (JSON only).

Cluster designation already rides the screen: ``gspec.cluster/2`` and the
rich ``gspec.locus_map/1`` carry the evaluated ``phenotypes`` map.  This
reader turns those documents into the same :class:`TypingFileResult` shape
the gene path renders, so ``gspec typing`` can produce the concise
designation table without re-running minimap2 or loading the database.

TSV cluster results are intentionally not parsed here: the per-scheme table
does not carry the database name or scheme metadata; JSON does.
"""

from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field, ValidationError

from gspec.errors import InputError
from gspec.formats.typing_result import TypingFileResult
from gspec.typing_results import SchemeCall

_CLUSTER_SCHEMAS = {"gspec.cluster/2", "gspec.locus_map/1"}


class _ClusterParams(BaseModel, frozen=True):
    db: str


class _ClusterFile(BaseModel, frozen=True):
    file: str
    phenotypes: dict[str, SchemeCall] = Field(default_factory=dict)


class _ClusterDocument(BaseModel, frozen=True):
    model_config = ConfigDict(populate_by_name=True)

    schema_name: str = Field(alias="schema")
    params: _ClusterParams
    files: list[_ClusterFile]


def _read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise InputError(
            f"cluster screen result not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        ) from exc


def _looks_like_json(paths: list[Path], stdin_text: str | None) -> bool:
    if stdin_text is not None:
        return stdin_text.lstrip().startswith("{")
    return bool(paths) and all(path.suffix.lower() == ".json" for path in paths)


def _parse_document(label: str, text: str) -> tuple[str, list[TypingFileResult]]:
    try:
        document = _ClusterDocument.model_validate_json(text)
    except ValidationError as exc:
        raise InputError(
            f"cluster screen result is not a valid gspec cluster/2 or locus_map/1 document:"
            f" {label}",
            code="TYPING_MALFORMED",
            context={"file": label},
        ) from exc
    if document.schema_name not in _CLUSTER_SCHEMAS:
        raise InputError(
            f"unsupported cluster result schema {document.schema_name!r}: {label}",
            code="TYPING_MALFORMED",
            context={"file": label, "schema": document.schema_name},
        )
    if not document.files:
        raise InputError(
            f"cluster screen result has no files to type: {label}",
            code="TYPING_NO_DATA",
            context={"file": label},
        )
    return document.params.db, [
        TypingFileResult(file=entry.file, genes=(), phenotypes=entry.phenotypes)
        for entry in document.files
    ]


def read_cluster_typing_input(
    paths: list[Path], stdin_text: str | None = None
) -> tuple[str, list[TypingFileResult]] | None:
    """Read cluster/2 or locus_map JSON, or return None for table input.

    ``None`` means the caller should fall back to the gene-table parser; a
    JSON-looking input with the wrong schema is a typed error instead.
    """
    if not _looks_like_json(paths, stdin_text):
        return None
    documents: list[tuple[str, str]] = [(str(path), _read_text(path)) for path in paths]
    if stdin_text is not None:
        documents.append(("-", stdin_text))
    databases: list[str] = []
    results: list[TypingFileResult] = []
    for label, text in documents:
        db, parsed = _parse_document(label, text)
        if db not in databases:
            databases.append(db)
        results.extend(parsed)
    if len(databases) > 1:
        raise InputError(
            "cluster screen results mix databases; type each database's results separately:"
            f" {', '.join(databases)}",
            code="DATABASE_MISMATCH",
            context={"databases": ", ".join(databases)},
        )
    return databases[0], results
