#!/usr/bin/env python3
"""Cluster-typing calibration harness — a DEVELOPER tool, not a CLI surface.

Runs the cluster screen for every labeled sample, scores the database's
gapit.typing/1 rules, and prints three reports for tuning cutoff /
ambiguity_margin / rule weights (the v1 loop the future learned_linear
training flow will consume):

  1. per-expected-phenotype score distribution (n/min/median/max over the
     expected phenotype's best rule score),
  2. the agreement matrix called x expected,
  3. the divergence list (sample, expected, called, scores).

Usage:
    python scripts/cluster_calibration.py DB --datadir DIR LABELS.tsv

LABELS.tsv is two columns, sample_path<TAB>expected_phenotype, one row per
sample (a leading '#' starts a comment; blank lines are skipped). The
database must be a typed cluster database (typing.json installed via
`gapit db build --typing`). A sample that fails to screen aborts with its
typed error; otherwise the exit code is 0 however badly the model scores.
"""

import argparse
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path

from gapit import config
from gapit.cluster import ClusterParams, load_features, load_typing, screen_cluster_file
from gapit.screening import find_database
from gapit.typing_engine import decide_typing, score_rules

NO_CALL = "(no call)"
AMBIGUOUS = "-"


@dataclass(frozen=True, slots=True)
class SampleResult:
    """One screened sample: the label, the call, and both scores."""

    sample: str
    expected: str
    called: str
    called_score: float
    expected_score: float | None
    confidence: str


def calibrate(
    db_name: str, datadir: Path | None, labels: list[tuple[str, str]], **thresholds: float
) -> list[SampleResult]:
    """Screen every labeled sample against the typed cluster db and reduce
    each to a calibration row (the importable seam the self-test drives)."""
    database = find_database(config.resolve_datadir(datadir), db_name)
    typing_document = load_typing(database)
    if typing_document is None:
        raise SystemExit(
            f"error: database {db_name} carries no typing.json (calibration needs a typed db)"
        )
    params = ClusterParams(db=db_name, **thresholds)
    results: list[SampleResult] = []
    for sample, expected in labels:
        report = screen_cluster_file(Path(sample), database, load_features(database), params)
        scored = score_rules(report, typing_document)
        phenotype, detail = decide_typing(scored, typing_document)
        if report.best is None:
            called, called_score, confidence = NO_CALL, 0.0, "no-locus"
        else:
            called = phenotype if phenotype is not None else AMBIGUOUS
            called_score, confidence = detail.score, detail.confidence
        expected_score = max(
            (rule.score for rule in scored if rule.phenotype == expected), default=None
        )
        results.append(
            SampleResult(
                sample=sample,
                expected=expected,
                called=called,
                called_score=called_score,
                expected_score=expected_score,
                confidence=confidence,
            )
        )
    return results


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="calibrate a typed cluster database against labels"
    )
    parser.add_argument("db", help="cluster database name under the datadir")
    parser.add_argument(
        "--datadir", type=Path, default=None, help="datadir (default: gapit resolution)"
    )
    parser.add_argument("labels", type=Path, help="TSV: sample_path<TAB>expected_phenotype")
    parser.add_argument("--min-gene-cov", type=float, default=90.0)
    parser.add_argument("--min-gene-id", type=float, default=90.0)
    parser.add_argument("--min-cluster-cov", type=float, default=96.0)
    args = parser.parse_args(argv)

    labels = read_labels(args.labels)
    if not labels:
        print("error: labels file has no rows", file=sys.stderr)
        return 5
    results = calibrate(
        args.db,
        args.datadir,
        labels,
        min_gene_cov=args.min_gene_cov,
        min_gene_id=args.min_gene_id,
        min_cluster_cov=args.min_cluster_cov,
    )
    print_distribution(results)
    print_matrix(results)
    print_divergences(results)
    return 0


def read_labels(path: Path) -> list[tuple[str, str]]:
    """sample<TAB>expected rows; '#' comments and blanks skipped."""
    rows: list[tuple[str, str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        sample, tab, expected = line.partition("\t")
        if tab and expected.strip():
            rows.append((sample.strip(), expected.strip()))
    return rows


def print_distribution(results: list[SampleResult]) -> None:
    """(a) per-expected-phenotype distribution of the expected rule's score."""
    by_expected: dict[str, list[float]] = {}
    for result in results:
        if result.expected_score is not None:
            by_expected.setdefault(result.expected, []).append(result.expected_score)
    print("== score distribution per expected phenotype (expected rule's score) ==")
    print("expected\tn\tmin\tmedian\tmax")
    for expected in sorted(by_expected):
        scores = by_expected[expected]
        print(
            f"{expected}\t{len(scores)}\t{min(scores):.4f}"
            f"\t{statistics.median(scores):.4f}\t{max(scores):.4f}"
        )
    unmatched = sorted({r.expected for r in results if r.expected_score is None})
    if unmatched:
        print(f"(no rule for expected phenotype: {', '.join(unmatched)})")


def print_matrix(results: list[SampleResult]) -> None:
    """(b) the agreement matrix: rows expected, columns called."""
    called_values = sorted({result.called for result in results})
    cells: dict[tuple[str, str], int] = {}
    for result in results:
        cells[(result.expected, result.called)] = cells.get((result.expected, result.called), 0) + 1
    print("\n== agreement matrix (rows: expected, columns: called) ==")
    print("expected\t" + "\t".join(called_values))
    for expected in sorted({result.expected for result in results}):
        row = "\t".join(str(cells.get((expected, called), 0)) for called in called_values)
        print(f"{expected}\t{row}")


def print_divergences(results: list[SampleResult]) -> None:
    """(c) every sample whose call disagrees with its label."""
    divergent = [result for result in results if result.called != result.expected]
    print(f"\n== divergences ({len(divergent)} of {len(results)}) ==")
    print("sample\texpected\tcalled\tcalled_score\texpected_score\tconfidence")
    for result in divergent:
        expected_score = "-" if result.expected_score is None else f"{result.expected_score:.4f}"
        print(
            f"{result.sample}\t{result.expected}\t{result.called}\t"
            f"{result.called_score:.4f}\t{expected_score}\t{result.confidence}"
        )


if __name__ == "__main__":
    sys.exit(main())
