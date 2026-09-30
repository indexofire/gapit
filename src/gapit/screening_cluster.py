"""The cluster-screening use-case and db-kind guards (stage 2/3).

Split from screening.py so each use-case module stays under the 250 pure-LOC
ceiling: this module owns the cluster branch — the shared CLI+MCP kind
guards, the gated ClusterParams constructor, and the orchestration that
screens each file and renders once. Stage 3 adds the typing wiring: a
cluster db with a ``typing.json`` is reference-validated up front
(TYPING_UNKNOWN_GENE) and every report is annotated with the phenotype call
before rendering; untyped databases render the stage-2 output unchanged.
"""

from datetime import UTC, datetime
from pathlib import Path

import typer

from gapit import db
from gapit.cluster import (
    ClusterParams,
    ClusterReport,
    load_features,
    load_typing,
    screen_cluster_file,
)
from gapit.db import Database
from gapit.engines import AlignerEnum, OutputFormat
from gapit.errors import usage_fail
from gapit.formats.cluster import format_cluster_tsv, render_cluster_json
from gapit.formats.cluster_md import render_cluster_md
from gapit.typing_engine import evaluate_typing
from gapit.typing_models import validate_references

_CLUSTER_DEFAULTS = (90.0, 90.0, 96.0)


def reject_gene_engine_flags(
    minid: float, mincov: float, merge_fragments: bool, jobs: int, aligner: AlignerEnum | None
) -> None:
    """The shared db-kind guard (CLI and MCP): any non-default gene-engine
    flag against a cluster database is a usage error naming the cluster
    flags — symmetric with :func:`reject_cluster_engine_flags`."""
    if minid != 80.0 or mincov != 80.0:
        usage_fail(
            "--minid/--mincov apply to gene databases only;"
            " cluster screening uses --min-gene-cov/--min-gene-id"
        )
    if merge_fragments:
        usage_fail("--merge-fragments is not available with a cluster database")
    if jobs != 1:
        usage_fail("--jobs is not available with a cluster database")
    if aligner is not None:
        usage_fail(
            "--aligner is not available with a cluster database"
            " (the cluster engine always uses minimap2)"
        )


def reject_cluster_engine_flags(
    min_gene_cov: float, min_gene_id: float, min_cluster_cov: float
) -> None:
    """Any non-default cluster flag on a gene database (or in reads mode) is
    a usage error — flags are never silently ignored."""
    if (min_gene_cov, min_gene_id, min_cluster_cov) != _CLUSTER_DEFAULTS:
        usage_fail("--min-gene-cov/--min-gene-id/--min-cluster-cov apply to cluster databases only")


def resolve_cluster_params(
    db: str, min_gene_cov: float, min_gene_id: float, min_cluster_cov: float, threads: int
) -> ClusterParams:
    """Usage-gated constructor: each cluster threshold is range-checked with
    its CLI flag name before the model is built (threads is validated
    upstream together with the gene-path checks)."""
    for flag, value in (
        ("--min-gene-cov", min_gene_cov),
        ("--min-gene-id", min_gene_id),
        ("--min-cluster-cov", min_cluster_cov),
    ):
        if not 0.0 <= value <= 100.0:
            usage_fail(f"{flag} must be in [0, 100]: got {value}")
    return ClusterParams(
        db=db,
        min_gene_cov=min_gene_cov,
        min_gene_id=min_gene_id,
        min_cluster_cov=min_cluster_cov,
        threads=threads,
    )


def _typed_reports(
    files: list[Path], database: Database, params: ClusterParams, *, quiet: bool, debug: bool
) -> tuple[list[ClusterReport], bool]:
    """Screen every file in order; annotate each report with the phenotype
    call when the db carries a typing.json (reference-validated once, up
    front). Returns the reports and whether the run is typed."""
    features = load_features(database)
    typing_document = load_typing(database)
    if typing_document is not None:
        validate_references(typing_document, features)
    reports: list[ClusterReport] = []
    for path in files:
        if not quiet:
            typer.echo(f"Processing: {path}", err=True)
        report = screen_cluster_file(path, database, features, params, debug=debug)
        if typing_document is not None:
            report = evaluate_typing(report, typing_document)
        if not quiet:
            called = report.best.locus if report.best is not None else "none"
            typer.echo(f"Best locus in {path}: {called}", err=True)
        reports.append(report)
    return reports, typing_document is not None


def run_cluster_screen(
    files: list[Path],
    database: db.Database,
    params: ClusterParams,
    *,
    output_format: OutputFormat,
    noheader: bool,
    nopath: bool,
    quiet: bool,
    debug: bool = False,
) -> str:
    """Cluster-engine orchestration beside the gene path: screen each input
    file in order, then render once (json/md/tsv via formats.cluster*)."""
    reports, typed = _typed_reports(files, database, params, quiet=quiet, debug=debug)
    now = datetime.now(UTC)
    if output_format is OutputFormat.json:
        return render_cluster_json(reports, params, now=now, typed=typed)
    if output_format is OutputFormat.md:
        return render_cluster_md(reports, params, now=now, typed=typed)
    return format_cluster_tsv(
        reports,
        csv=output_format is OutputFormat.csv,
        noheader=noheader,
        nopath=nopath,
        typed=typed,
    )
