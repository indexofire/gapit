"""Typed errors, shared raise-helpers, and the JSON error envelope
(gapit.error/1).

Exit-code contract (AGENTS.md §5): 2 usage, 3 missing dependency, 4 db error,
5 input error, 1 unexpected. Failures render as a one-line JSON envelope on
stderr via ``render_error``.
"""

from pathlib import Path
from typing import Literal, NoReturn, cast

from pydantic import BaseModel, ConfigDict, Field, JsonValue

# Public error-code vocabulary. Keep in sync with every raise site; the test
# suite AST-scans the source so a new code cannot be added without registering
# it here.
ERROR_CODES: frozenset[str] = frozenset(
    {
        "ARCHIVE_INVALID",
        "ARCHIVE_TOO_LARGE",
        "BLASTDBCMD_PARSE_FAILED",
        "BLAST_FAILED",
        "BLAST_PARSE_FAILED",
        "BLAST_VERSION_CHECK_FAILED",
        "BLAST_VERSION_TOO_OLD",
        "BUILD_INVALID",
        "BUILD_SELF_CHECK_FAILED",
        "BUNDLED_MALFORMED",
        "BUNDLED_NOT_FOUND",
        "CHECKSUM_MISMATCH",
        "CS_MISSING",
        "CS_PARSE_FAILED",
        "DATABASE_MISMATCH",
        "DATABASE_NOT_FOUND",
        "DATABASE_NOT_INDEXED",
        "DATADIR_CREATE_FAILED",
        "DATADIR_EMPTY",
        "DATADIR_NOT_FOUND",
        "DB_ALREADY_EXISTS",
        "DB_INCOMPLETE",
        "DBINFO_MALFORMED",
        "DB_LOCK_FAILED",
        "DOWNLOAD_FAILED",
        "DOWNLOAD_TOO_LARGE",
        "FEATURES_MALFORMED",
        "FEATURES_MISSING",
        "FILE_IO_ERROR",
        "FLOORS_MALFORMED",
        "FLOORS_UNKNOWN_GENE",
        "GENE_OUT_OF_BOUNDS",
        "GENE_WITHOUT_NAME",
        "GFF_WITHOUT_FASTA",
        "GSEEK_FAILED",
        "GSEEK_RECORD_INVALID",
        "HEADER_MALFORMED",
        "INPUT_NOT_FOUND",
        "INVALID_ENGINE",
        "INVALID_FASTA",
        "INVALID_GENBANK",
        "INVALID_GFF",
        "INVALID_INPUT",
        "INVALID_READS_FORMAT",
        "LOCUS_WITHOUT_GENES",
        "MAKEBLASTDB_FAILED",
        "MANIFEST_MALFORMED",
        "METADATA_MALFORMED",
        "MINIMAP2_FAILED",
        "MISSING_DEPENDENCY",
        "PAF_PARSE_FAILED",
        "PROVIDER_EMPTY",
        "PROVIDER_INVALID",
        "PROVIDER_MALFORMED",
        "RECORDS_MALFORMED",
        "SCREEN_TABLE_MALFORMED",
        "SUMMARY_MALFORMED",
        "TYPING_CLUSTER_DB",
        "TYPING_MALFORMED",
        "TYPING_NO_DATA",
        "TYPING_NO_SCHEME",
        "TYPING_UNKNOWN_GENE",
        "UNEXPECTED",
        "UNSUPPORTED_LOCATION",
        "USAGE_ERROR",
    }
)

# list is invariant; values are plain strings at this JSON-schema boundary.
ERROR_CODE_ENUM = cast(list[JsonValue], sorted(ERROR_CODES))


class GapitError(Exception):
    """Base class for all typed gapit errors."""

    code: str
    exit_code: int = 1

    def __init__(self, message: str, *, code: str, context: dict[str, str] | None = None) -> None:
        super().__init__(message)
        self.code = code
        self.context: dict[str, str] = context if context is not None else {}


class UsageError(GapitError):
    """Invalid command-line usage."""

    exit_code = 2


class DependencyError(GapitError):
    """An external binary (BLAST+, minimap2) is missing from PATH."""

    exit_code = 3


class DatabaseError(GapitError):
    """A datadir/database is missing, unindexed, or failed to build."""

    exit_code = 4


class InputError(GapitError):
    """User-supplied input is malformed."""

    exit_code = 5


def usage_fail(message: str) -> NoReturn:
    """Raise a usage error (gapit.error/1 envelope, exit 2)."""
    raise UsageError(message, code="USAGE_ERROR")


def ensure_input_file(path: Path, what: str = "input file") -> None:
    """Raise INPUT_NOT_FOUND for a missing/unreadable input path; ``what``
    names the kind in the message (e.g. "reads file")."""
    if not path.is_file():
        raise InputError(
            f"{what} not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        )


class ErrorEnvelope(BaseModel, frozen=True):
    """The gapit.error/1 stderr envelope."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.error/1"] = Field(default="gapit.error/1", alias="schema")
    code: str = Field(json_schema_extra={"enum": ERROR_CODE_ENUM})
    message: str
    context: dict[str, str]


def render_error(exc: BaseException) -> str:
    """One-line compact JSON envelope for any failure (non-GapitError becomes
    UNEXPECTED, exit 1)."""
    if isinstance(exc, GapitError):
        envelope = ErrorEnvelope(code=exc.code, message=str(exc), context=exc.context)
    else:
        envelope = ErrorEnvelope(
            code="UNEXPECTED", message=f"{type(exc).__name__}: {exc}", context={}
        )
    return envelope.model_dump_json(by_alias=True)
