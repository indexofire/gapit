"""Reads-mode CLI argument rules for `gapit screen`.

Split from cmd_screen.py at the 250 pure-LOC ceiling (the engines.py
precedent): the typer surface stays there; this module owns the
`--r1`/`--r2` comma-list decoding and the reads-mode flag rejects whose
exact messages are frozen by tests. The engine use-cases live in
screening_reads.py.
"""

from pathlib import Path

from gapit.errors import usage_fail


def split_read_list(raw: str, flag: str) -> list[Path]:
    """Split a comma-separated --r1/--r2 value into paths; empty elements
    are usage errors (an empty string would silently become the cwd)."""
    parts = [part.strip() for part in raw.split(",")]
    if any(not part for part in parts):
        usage_fail(f"{flag} contains an empty element: {raw!r}")
    return [Path(part) for part in parts]


def reject_reads_mode_flags(
    fofn: Path | None, noheader: bool, nopath: bool, jobs: int, merge_fragments: bool
) -> None:
    """--fofn/--noheader/--nopath/--jobs/--merge-fragments are blastn-contig
    flags; reads mode rejects each instead of silently ignoring it."""
    if fofn is not None:
        usage_fail("--fofn is not available in reads mode")
    if noheader:
        usage_fail("--noheader is not available in reads mode")
    if nopath:
        usage_fail("--nopath is not available in reads mode")
    if jobs != 1:
        usage_fail("--jobs is not available in reads mode")
    if merge_fragments:
        usage_fail("--merge-fragments is not available in reads mode")
