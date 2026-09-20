"""Shared external-tool plumbing: the argv subprocess runner and stderr notes."""

import subprocess
import sys

from gapit.errors import DependencyError


def run_tool(argv: list[str]) -> subprocess.CompletedProcess[str]:
    """Run an external tool with an argv list (never a shell)."""
    try:
        return subprocess.run(argv, check=False, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise DependencyError(
            f"required binary not found on PATH: {argv[0]}",
            code="MISSING_DEPENDENCY",
            context={"binary": argv[0]},
        ) from exc


def note(quiet: bool, message: str) -> None:
    """Per-step progress on stderr when quiet is disabled (stdout stays pure)."""
    if not quiet:
        print(f"gapit: {message}", file=sys.stderr)
