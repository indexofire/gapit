"""Positional-FASTQ sample auto-pairing for `gapit screen *.fastq.gz`.

The wildcard workflow (rightsholder request): a shell glob hands many gz
FASTQ files to ``screen`` and gapit infers samples from filename
conventions. Pure filename arithmetic — no I/O, no sniffing (the reads
engine owns content detection via ``gapit.reads.detect_read_kind``).

Conventions, matched case-insensitively, longest suffix first:
``_R1_001``/``_R2_001`` (bcl2fastq; a ``_L001`` lane tag stays part of the
sample key), ``_R1``/``_R2``, ``_1``/``_2``, ``.1``/``.2``. Extensions
``.fastq``/``.fq`` optionally under ``.gz``/``.bz2`` strip before marker
detection. A file with no detectable marker, or whose mate is missing,
screens single-end with a warning — never an error. Duplicate identical
basenames (two lane directories) merge as extra lanes of one sample and
union at the sample level.
"""

from collections.abc import Callable
from pathlib import Path

from pydantic import BaseModel

_READ_EXTENSIONS = (".fastq", ".fq")
_COMPRESSION_EXTENSIONS = (".gz", ".bz2")

# Longest suffix first: _R1_001 must win over _R1 before _1. Side 0 marks a
# marker-less file (its stem is the sample key; it never pairs).
_MATE_MARKERS: tuple[tuple[str, int], ...] = (
    ("_r1_001", 1),
    ("_r2_001", 2),
    ("_r1", 1),
    ("_r2", 2),
    ("_1", 1),
    ("_2", 2),
    (".1", 1),
    (".2", 2),
)

_UNMARKED = 0


class SampleLanes(BaseModel, frozen=True):
    """One inferred sample: its key plus its (r1, r2) lanes (r2 None =
    single-end). Lanes union into one files[] entry at the document level."""

    sample: str
    lanes: tuple[tuple[Path, Path | None], ...]


def _strip_read_extensions(name: str) -> str:
    """`sample_1.fastq.gz` -> `sample_1` (one compression suffix under one
    read extension, case-insensitive; unknown extensions pass through)."""
    stem = name
    lowered = stem.lower()
    for compression in _COMPRESSION_EXTENSIONS:
        if lowered.endswith(compression):
            stem = stem[: -len(compression)]
            break
    lowered = stem.lower()
    for extension in _READ_EXTENSIONS:
        if lowered.endswith(extension):
            return stem[: -len(extension)]
    return stem


def _split_mate_marker(stem: str) -> tuple[str, int] | None:
    """`sample_R1_001` -> (`sample`, 1); None when no marker matches."""
    lowered = stem.lower()
    for marker, side in _MATE_MARKERS:
        if lowered.endswith(marker):
            return stem[: -len(marker)], side
    return None


def pair_samples(paths: list[Path], *, warn: Callable[[str], None]) -> list[SampleLanes]:
    """Infer samples from reads filenames. Within a sample, R1-marked and
    R2-marked files pair in sorted (name, path) order; surplus files of
    either side, and marker-less files, become single-end lanes, each with a
    warning. Samples sort lexicographically, lanes by original name."""
    grouped: dict[str, dict[int, list[Path]]] = {}
    for path in paths:
        stem = _strip_read_extensions(path.name)
        split = _split_mate_marker(stem)
        if split is None:
            grouped.setdefault(stem, {}).setdefault(_UNMARKED, []).append(path)
            continue
        sample, side = split
        grouped.setdefault(sample, {}).setdefault(side, []).append(path)
    samples: list[SampleLanes] = []
    for sample in sorted(grouped):
        sides = grouped[sample]
        r1s = sorted(sides.get(1, []), key=lambda path: (path.name, str(path)))
        r2s = sorted(sides.get(2, []), key=lambda path: (path.name, str(path)))
        unmarked = sorted(sides.get(_UNMARKED, []), key=lambda path: (path.name, str(path)))
        lanes: list[tuple[Path, Path | None]] = []
        for index in range(max(len(r1s), len(r2s))):
            r1 = r1s[index] if index < len(r1s) else None
            r2 = r2s[index] if index < len(r2s) else None
            if r1 is not None and r2 is not None:
                lanes.append((r1, r2))
                continue
            lone = r1 if r1 is not None else r2
            if lone is not None:  # narrowed for the type checker
                warn(f"no mate found for {lone} — screening single-end")
                lanes.append((lone, None))
        for path in unmarked:
            warn(f"no mate found for {path} — screening single-end")
            lanes.append((path, None))
        samples.append(SampleLanes(sample=sample, lanes=tuple(lanes)))
    return samples
