"""Typed errors carrying machine-readable codes and documented exit codes.

Exit-code contract (AGENTS.md §5): 2 usage, 3 missing dependency, 4 db error,
5 input error, 1 unexpected. The JSON error envelope lands in Phase 4; the CLI
currently renders ``ERROR: <message>`` on stderr.
"""


class GaitaError(Exception):
    """Base class for all typed gaita errors."""

    code: str
    exit_code: int = 1

    def __init__(self, message: str, *, code: str) -> None:
        super().__init__(message)
        self.code = code


class DependencyError(GaitaError):
    """An external binary (BLAST+, any2fasta) is missing from PATH."""

    exit_code = 3


class DatabaseError(GaitaError):
    """A datadir/database is missing, unindexed, or failed to build."""

    exit_code = 4


class InputError(GaitaError):
    """User-supplied input is malformed."""

    exit_code = 5
