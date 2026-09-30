"""Target-space cs math for the cluster screening engine (stage 2 core).

minimap2's short-form ``cs`` tag describes one alignment as a walk of ops;
this module turns that walk into half-open target intervals (``CsRun``) and
reduces runs across records to the (covered, matched) base counts behind
per-gene and per-locus coverage/identity. Pure functions only — no I/O, no
minimap2; the engine (gapit.cluster) owns invocation and verdict assembly.
"""

from collections.abc import Iterable
from typing import Literal

from gapit.errors import GapitError
from gapit.paf import union_length

# One cs op projected onto the target: half-open [start, end) with `matches`
# matching bases inside it (substitution and deletion runs carry 0).
CsRun = tuple[int, int, int]

Verdict = Literal["present", "partial", "absent"]

# Verdict floor shared by every caller (documented constant): a gene covered
# at >= 50% is "partial" once it fails either present threshold, below it
# "absent".
PARTIAL_MIN_COV = 50.0


def parse_cs(cs: str, tstart: int) -> tuple[CsRun, ...]:
    """Parse a short-form cs string starting at PAF ``tstart`` into target
    runs. ``:n`` (n matching bases), ``*ab`` (substitution: 1 target base, 0
    matches), ``-ss`` (deletion: len(ss) target bases, 0 matches) consume
    target columns; ``+s`` (insertion) consumes none. Any other op shape is a
    typed CS_PARSE_FAILED error (the PAF boundary's parse-don't-validate)."""
    runs: list[CsRun] = []
    position = tstart
    index = 0
    while index < len(cs):
        op = cs[index]
        index += 1
        if op == ":":
            digits_start = index
            while index < len(cs) and cs[index].isdigit():
                index += 1
            if digits_start == index:
                raise GapitError(
                    f"cs match op without an integer length: {cs!r}",
                    code="CS_PARSE_FAILED",
                )
            length = int(cs[digits_start:index])
            runs.append((position, position + length, length))
            position += length
        elif op == "*":
            if index + 2 > len(cs):
                raise GapitError(
                    f"cs substitution op without two bases: {cs!r}", code="CS_PARSE_FAILED"
                )
            runs.append((position, position + 1, 0))
            position += 1
            index += 2
        elif op == "+":
            bases_start = index
            while index < len(cs) and cs[index] not in ":*+-":
                index += 1
            if bases_start == index:
                raise GapitError(f"cs insertion op without bases: {cs!r}", code="CS_PARSE_FAILED")
            # insertion bases: query-only, target position unchanged
        elif op == "-":
            bases_start = index
            while index < len(cs) and cs[index] not in ":*+-":
                index += 1
            if bases_start == index:
                raise GapitError(f"cs deletion op without bases: {cs!r}", code="CS_PARSE_FAILED")
            length = index - bases_start
            runs.append((position, position + length, 0))
            position += length
        else:
            raise GapitError(
                f"unknown cs op {op!r} in {cs!r}",
                code="CS_PARSE_FAILED",
            )
    return tuple(runs)


def intersect_coverage(records: Iterable[Iterable[CsRun]], lo: int, hi: int) -> tuple[int, int]:
    """(covered, matched) target bases of the half-open [lo, hi) across every
    record's runs: covered is the union of all run slices (de-duplicated per
    base, so cross-contig fragments union and overlaps never double-count);
    matched is the union of the matching slices — a base counts as matched
    when ANY record matching run covers it. Coverage denominators come from
    the caller; identity is matched/covered (0 when covered is 0)."""
    match_intervals: list[tuple[int, int]] = []
    all_intervals: list[tuple[int, int]] = []
    for runs in records:
        for start, end, matches in runs:
            if end <= lo or start >= hi:
                continue
            clipped = (max(start, lo), min(end, hi))
            all_intervals.append(clipped)
            if matches > 0:
                match_intervals.append(clipped)
    return union_length(all_intervals), union_length(match_intervals)


def gene_verdict(
    coverage_pct: float,
    identity_pct: float,
    *,
    min_gene_cov: float,
    min_gene_id: float,
) -> Verdict:
    """present (both thresholds met, boundary-inclusive) / partial (covered
    at >= PARTIAL_MIN_COV but failing a present threshold) / absent."""
    if coverage_pct >= min_gene_cov and identity_pct >= min_gene_id:
        return "present"
    if coverage_pct >= PARTIAL_MIN_COV:
        return "partial"
    return "absent"


def locus_sort_key(
    locus_id: str, coverage_pct: float, identity_pct: float, covered_bp: int
) -> tuple[float, float, int, str]:
    """Rank key: coverage desc, identity desc, covered bp desc (short-
    reference bias guard), locus id asc (determinism). Sort ascending."""
    return (-coverage_pct, -identity_pct, -covered_bp, locus_id)
