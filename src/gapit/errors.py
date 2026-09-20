"""Typed errors and the JSON error envelope (gapit.error/1).

Exit-code contract (AGENTS.md §5): 2 usage, 3 missing dependency, 4 db error,
5 input error, 1 unexpected. Failures render as a one-line JSON envelope on
stderr via ``render_error``.
"""

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field


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


class ErrorEnvelope(BaseModel, frozen=True):
    """The gapit.error/1 stderr envelope."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.error/1"] = Field(default="gapit.error/1", alias="schema")
    code: str
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
