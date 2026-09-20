"""Property tests: interval-union coverage equals the per-base oracle.

The pre-optimization aggregate_coverage allocated one Python int per target
base and incremented it per alignment position. These tests keep that
implementation as a local oracle and assert the interval-union rewrite is
exactly equivalent (identical integers in, identical floats out) across
edge cases, seeded fuzz, and a 50k-row smoke benchmark.
"""

import random
import time
from collections.abc import Mapping, Sequence
from typing import Literal

from gapit.dbcodec import decode_seqid
from gapit.paf import PafRecord, parse_paf_row, union_length
from gapit.reads import GeneCoverage, aggregate_coverage

Interval = tuple[int, int]

# ---------------------------------------------------------------------------
# Oracles: the OLD (pre-rewrite) implementations, copied verbatim.
# ---------------------------------------------------------------------------


def per_base_covered(intervals: Sequence[Interval], tlen: int) -> int:
    """Oracle: covered bases by counting positions with depth > 0."""
    depth = [0] * tlen
    for start, end in intervals:
        for position in range(start, end):
            depth[position] += 1
    return sum(1 for value in depth if value > 0)


def aggregate_coverage_oracle(
    rows: list[PafRecord],
    *,
    default_db: str,
    min_breadth: float,
    products: Mapping[str, str] | None = None,
) -> list[GeneCoverage]:
    """Oracle: the pre-rewrite per-base aggregate_coverage, copied, plus the
    reads/2 alen-weighted identity mean restated naively."""
    depths: dict[str, list[int]] = {}
    qnames: dict[str, set[str]] = {}
    weights: dict[str, list[tuple[float, int]]] = {}
    for row in rows:
        if not row.is_primary:
            continue
        if row.tname not in depths:
            depths[row.tname] = [0] * row.tlen
            qnames[row.tname] = set()
            weights[row.tname] = []
        depth = depths[row.tname]
        for position in range(row.tstart, row.tend):
            depth[position] += 1
        qnames[row.tname].add(row.qname)
        identity = (
            100.0 * (row.alen - row.nm) / row.alen if row.nm is not None and row.alen > 0 else 100.0
        )
        weights[row.tname].append((identity, row.alen))
    descriptions = products or {}
    coverages: list[GeneCoverage] = []
    for tname, depth in depths.items():
        tlen = len(depth)
        covered = sum(1 for value in depth if value > 0)
        header = decode_seqid(tname, default_db)
        pairs = weights[tname]
        weight_sum = sum(weight for _, weight in pairs)
        coverages.append(
            GeneCoverage(
                database=header.database,
                gene=header.gene,
                accession=header.accession,
                function=header.function,
                product=descriptions.get(tname, ""),
                tlen=tlen,
                breadth_pct=100.0 * covered / tlen,
                mean_depth=sum(depth) / tlen,
                reads_mapped=len(qnames[tname]),
                present=100.0 * covered / tlen >= min_breadth,
                mean_identity_pct=(
                    sum(identity * weight for identity, weight in pairs) / weight_sum
                    if weight_sum
                    else 0.0
                ),
            )
        )
    return sorted(coverages, key=lambda entry: (-entry.breadth_pct, entry.gene))


def paf(
    qname: str,
    tname: str,
    tlen: int,
    tstart: int,
    tend: int,
    *,
    is_primary: bool = True,
    nm: int | None = None,
) -> PafRecord:
    """A PAF row spanning [tstart, tend) on a target of length tlen."""
    return PafRecord(
        qname=qname,
        qlen=tend - tstart,
        qstart=0,
        qend=tend - tstart,
        strand="+",
        tname=tname,
        tlen=tlen,
        tstart=tstart,
        tend=tend,
        nmatch=tend - tstart,
        alen=tend - tstart,
        mapq=60,
        is_primary=is_primary,
        nm=nm,
    )


# ---------------------------------------------------------------------------
# union_length: fixed edge cases + seeded fuzz vs per-base counting.
# ---------------------------------------------------------------------------


def test_union_length_edge_cases_match_per_base() -> None:
    """Given hand-built interval sets (empty, single, nested, adjacent,
    disjoint, zero-length, unsorted, tlen-filling), When union_length runs,
    Then every result equals the per-base count oracle."""
    cases: list[tuple[int, list[Interval]]] = [
        (100, []),
        (100, [(0, 0)]),
        (100, [(0, 100)]),
        (100, [(23, 24)]),
        (100, [(0, 50), (25, 75)]),
        (100, [(0, 50), (10, 20)]),
        (100, [(0, 50), (50, 100)]),
        (100, [(0, 10), (90, 100)]),
        (100, [(3, 3), (0, 10)]),
        (100, [(0, 10), (3, 3), (3, 3)]),
        (100, [(0, 10), (5, 15), (12, 20)]),
        (100, [(50, 60), (0, 10), (30, 40)]),
        (1, [(0, 1)]),
        (100, [(100, 100)]),
    ]
    for tlen, intervals in cases:
        assert union_length(intervals) == per_base_covered(intervals, tlen), (tlen, intervals)


def test_union_length_fuzz_matches_per_base() -> None:
    """Given 200 seeded random interval sets (varying density and overlap),
    When union_length runs, Then covered bases are identical to the per-base
    oracle in every case."""
    rng = random.Random(20260920)
    for _ in range(200):
        tlen = rng.randint(1, 500)
        intervals: list[Interval] = []
        for _ in range(rng.randint(0, 80)):
            start = rng.randint(0, tlen)
            end = rng.randint(start, tlen)
            intervals.append((start, end))
        assert union_length(intervals) == per_base_covered(intervals, tlen), intervals


# ---------------------------------------------------------------------------
# aggregate_coverage: seeded fuzz vs the full old implementation.
# ---------------------------------------------------------------------------

TNAME_A = "db~~~geneA~~~ACC~~~RES"
TNAME_B = "db~~~geneB~~~ACC~~~RES"
TNAME_POOL = [TNAME_A, TNAME_B, "db~~~geneC~~~ACC~~~RES", "plainid"]


def test_aggregate_fuzz_matches_oracle() -> None:
    """Given 200 seeded random alignment sets (multi-target, overlapping,
    repeated qnames, secondary rows, inconsistent tlen across rows), When
    aggregated by the rewritten aggregate_coverage, Then every GeneCoverage
    field is identical to the per-base oracle — locking the tlen-first-seen
    denominator, breadth/mean arithmetic, qname sets, and sort order."""
    rng = random.Random(81021)
    for _ in range(200):
        base_tlen = {name: rng.randint(1, 200) for name in TNAME_POOL}
        rows: list[PafRecord] = []
        for _ in range(rng.randint(0, 60)):
            tname = rng.choice(TNAME_POOL)
            limit = base_tlen[tname]
            # Rows may carry a tlen larger than the base; the FIRST row's
            # tlen is the denominator in both implementations, and intervals
            # never exceed the base so the oracle's array never overflows.
            tlen = rng.choice([limit, limit + rng.randint(0, 50)])
            start = rng.randint(0, limit)
            end = rng.randint(start, limit)
            rows.append(
                paf(
                    qname=f"r{rng.randint(1, 8)}",
                    tname=tname,
                    tlen=tlen,
                    tstart=start,
                    tend=end,
                    is_primary=rng.random() > 0.2,
                    nm=rng.choice([None, None, rng.randint(0, max(0, end - start))]),
                )
            )
        assert aggregate_coverage(
            rows, default_db="db", min_breadth=90.0
        ) == aggregate_coverage_oracle(rows, default_db="db", min_breadth=90.0), rows


def test_aggregate_products_flow_through_unchanged() -> None:
    """Given rows plus a products mapping, When aggregated, Then product
    descriptions match the oracle exactly (same lookup contract)."""
    rows = [paf("r1", TNAME_A, 100, 0, 60), paf("r2", TNAME_B, 80, 10, 30)]
    products = {TNAME_A: "gene A product", TNAME_B: "gene B product"}
    assert aggregate_coverage(
        rows, default_db="db", min_breadth=50.0, products=products
    ) == aggregate_coverage_oracle(rows, default_db="db", min_breadth=50.0, products=products)


# ---------------------------------------------------------------------------
# Large-scale smoke: parse + aggregate 50k rows through both paths.
# ---------------------------------------------------------------------------


def _synthesize_paf(rng: random.Random, n_rows: int, n_targets: int, tlen: int) -> list[str]:
    """Deterministic PAF text: one line per alignment over n_targets genes."""
    lines: list[str] = []
    for i in range(n_rows):
        gene = i % n_targets
        start = rng.randint(0, tlen - 2)
        end = rng.randint(start + 1, min(tlen, start + 300))
        tags: Literal["tp:A:P", "tp:A:S"] = "tp:A:S" if i % 17 == 0 else "tp:A:P"
        fields: list[str] = [
            f"read{i % 9000}",
            str(end - start),
            "0",
            str(end - start),
            "+",
            f"db~~~gene{gene:04d}~~~ACC~~~RES",
            str(tlen),
            str(start),
            str(end),
            str(end - start),
            str(end - start),
            "60",
            tags,
        ]
        lines.append("\t".join(fields))
    return lines


def test_fifty_k_row_smoke_matches_oracle_and_benchmarks() -> None:
    """Given ~50k synthesized PAF rows, When parsed and aggregated through
    both the rewrite and the per-base oracle, Then the gene lists are
    identical and the wall-time delta is reported (benchmark evidence)."""
    rng = random.Random(4242)
    lines = _synthesize_paf(rng, n_rows=50_000, n_targets=500, tlen=1000)
    rows = [parse_paf_row(line) for line in lines]

    new_start = time.perf_counter()
    new = aggregate_coverage(rows, default_db="db", min_breadth=90.0)
    new_elapsed = time.perf_counter() - new_start

    oracle_start = time.perf_counter()
    old = aggregate_coverage_oracle(rows, default_db="db", min_breadth=90.0)
    oracle_elapsed = time.perf_counter() - oracle_start

    assert new == old
    assert len(new) == 500
    print(f"\naggregate 50k rows: union-sweep {new_elapsed:.3f}s vs per-base {oracle_elapsed:.3f}s")
