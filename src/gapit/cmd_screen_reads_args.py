"""Reads-mode CLI argument rules for `gapit screen`.

Split from cmd_screen.py at the 250 pure-LOC ceiling (the engines.py
precedent): the typer surface stays there; this module owns the
`--r1`/`--r2` comma-list decoding, the positional-FASTQ classification that
routes an all-FASTQ wildcard into reads mode, and the reads-mode flag
rejects whose exact messages are frozen by tests. The engine use-cases live
in screening_reads.py and screening_reads_positional.py.
"""

from pathlib import Path

from gapit.errors import InputError, usage_fail
from gapit.reads import ReadFileKind, detect_read_kind

# What the reads engine can actually take: raw or gzipped .fastq/.fq. A
# .fastq.bz2 name is NOT classified from extension — the sniff below cannot
# see through bzip2, so the file keeps the blastn path it has today.
_READS_SUFFIXES = (".fastq.gz", ".fastq", ".fq.gz", ".fq", ".fastq.gzip", ".fq.gzip")
# Contig-side extensions seqconvert owns; never sniffed, so the positional
# FASTA/GBK/EMBL path stays byte-identical.
_ASSEMBLY_SUFFIXES = (
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
    ".fsa.gz",
    ".gbk.gz",
    ".gb.gz",
    ".embl.gz",
    ".emb.gz",
    ".fa.bz2",
    ".fasta.bz2",
    ".fna.bz2",
    ".fsa.bz2",
    ".gbk.bz2",
    ".gb.bz2",
    ".embl.bz2",
    ".emb.bz2",
    ".fa",
    ".fasta",
    ".fna",
    ".fsa",
    ".gbk",
    ".gb",
    ".embl",
    ".emb",
)


def reads_positional(path: Path) -> bool:
    """Whether one positional FILE should route to reads mode: a .fastq/.fq
    (± .gz) extension says yes outright; known contig extensions say no
    outright; anything else is ambiguous and gets one content sniff through
    detect_read_kind (reused, never duplicated). Unreadable or
    undecipherable content keeps the contig path, which owns its typed
    input errors."""
    lowered = path.name.lower()
    if lowered.endswith(_READS_SUFFIXES):
        return True
    if lowered.endswith(_ASSEMBLY_SUFFIXES):
        return False
    try:
        return detect_read_kind(path) is ReadFileKind.fastq
    except (InputError, OSError):
        return False


def reject_mixed_positionals(files: list[Path], reads_flags: list[bool]) -> None:
    """A glob mixing reads and contig inputs is a usage error naming the
    reads files (the ones that pulled toward the reads engine)."""
    if not reads_flags or all(reads_flags):
        return
    offenders = ", ".join(
        str(path) for path, is_reads in zip(files, reads_flags, strict=True) if is_reads
    )
    usage_fail(f"mixed assembly and reads inputs; screen them separately: {offenders}")


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
