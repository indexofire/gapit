"""Target-space cs math for the cluster screening engine (stage 2 core).

minimap2's short-form ``cs`` tag describes one alignment as a walk of ops;
this module turns that walk into half-open target intervals (``CsRun``) and
reduces runs across records to the (covered, matched) base counts behind
per-gene and per-locus coverage/identity. Pure functions only — no I/O, no
minimap2; the engine (gapit.cluster) owns invocation and verdict assembly.
"""

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

from gapit.errors import GapitError
from gapit.paf import union_length

# One cs op projected onto the target: half-open [start, end) with `matches`
# matching bases inside it (substitution and deletion runs carry 0).
CsRun = tuple[int, int, int]

AlignmentOp = Literal["match", "substitution", "insertion", "deletion"]

Verdict = Literal["present", "partial", "absent"]


@dataclass(frozen=True, slots=True)
class CsStep:
    """One cs op in alignment order, with query and target coordinates
    relative to the alignment start (both 0-based half-open). Insertions
    consume query only, deletions target only; match/substitution consume
    both. ``matches`` is the number of matching target bases (0 for
    substitution/deletion/insertion)."""

    op: AlignmentOp
    qstart: int
    qend: int
    tstart: int
    tend: int
    matches: int


# Verdict floor shared by every caller (documented constant): a gene covered
# at >= 50% is "partial" once it fails either present threshold, below it
# "absent".
PARTIAL_MIN_COV = 50.0


def parse_cs_steps(cs: str) -> tuple[CsStep, ...]:
    """Parse a short-form cs string into query/target walk steps.

    ``:n`` (n matching bases), ``*ab`` (substitution: 1 query + 1 target
    base), ``+s`` (insertion: len(s) query bases), ``-ss`` (deletion:
    len(ss) target bases). Coordinates are relative to the alignment start.
    Any other op shape is a typed CS_PARSE_FAILED error (the PAF boundary's
    parse-don't-validate)."""
    steps: list[CsStep] = []
    qpos = 0
    tpos = 0
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
            steps.append(CsStep("match", qpos, qpos + length, tpos, tpos + length, matches=length))
            qpos += length
            tpos += length
        elif op == "*":
            if index + 2 > len(cs):
                raise GapitError(
                    f"cs substitution op without two bases: {cs!r}", code="CS_PARSE_FAILED"
                )
            steps.append(CsStep("substitution", qpos, qpos + 1, tpos, tpos + 1, matches=0))
            qpos += 1
            tpos += 1
            index += 2
        elif op == "+":
            bases_start = index
            while index < len(cs) and cs[index] not in ":*+-":
                index += 1
            if bases_start == index:
                raise GapitError(f"cs insertion op without bases: {cs!r}", code="CS_PARSE_FAILED")
            length = index - bases_start
            steps.append(CsStep("insertion", qpos, qpos + length, tpos, tpos, matches=0))
            qpos += length
        elif op == "-":
            bases_start = index
            while index < len(cs) and cs[index] not in ":*+-":
                index += 1
            if bases_start == index:
                raise GapitError(f"cs deletion op without bases: {cs!r}", code="CS_PARSE_FAILED")
            length = index - bases_start
            steps.append(CsStep("deletion", qpos, qpos, tpos, tpos + length, matches=0))
            tpos += length
        else:
            raise GapitError(
                f"unknown cs op {op!r} in {cs!r}",
                code="CS_PARSE_FAILED",
            )
    return tuple(steps)


def parse_cs(cs: str, tstart: int) -> tuple[CsRun, ...]:
    """Target-space projection of :func:`parse_cs_steps` (the historical
    coverage contract): match/substitution/deletion runs at absolute
    ``tstart``; insertions consume no target columns and are absent."""
    return tuple(
        (tstart + step.tstart, tstart + step.tend, step.matches)
        for step in parse_cs_steps(cs)
        if step.op != "insertion"
    )


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
