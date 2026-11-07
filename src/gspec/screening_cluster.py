"""The cluster-screening use-case and db-kind guards (stage 2/3).

Split from screening.py so each use-case module stays under the 250 pure-LOC
ceiling: this module owns the cluster branch — the shared CLI+MCP kind
guards, the gated ClusterParams constructor, and the orchestration that
screens each file and renders once. Stage 3 adds the typing wiring: a
cluster db with a ``typing.json`` is reference-validated up front
(TYPING_UNKNOWN_GENE) and every report is annotated before rendering — a
single-scheme document annotates the /1 best call (the byte-frozen stage-3
output), a multi-scheme document calls every scheme into the /2 phenotypes
map; untyped databases render the stage-2 output unchanged.
"""

import os
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path

import click

from gspec import db
from gspec.cluster import (
    ClusterParams,
    ClusterReport,
    load_features,
    load_typing,
    screen_cluster_file,
)
from gspec.cluster_profile import load_cluster_profile
from gspec.db import Database
from gspec.engines import AlignerEnum, Emit, OutputFormat
from gspec.errors import usage_fail
from gspec.formats.cluster import cluster_tsv_file_chunk, cluster_tsv_preamble, render_cluster_json
from gspec.formats.cluster2 import (
    cluster2_tsv_file_chunk,
    cluster2_tsv_preamble,
    render_cluster2_json,
)
from gspec.formats.cluster_md import (
    cluster2_md_file_chunk,
    cluster2_md_head,
    cluster_md_file_chunk,
    cluster_md_head,
)
from gspec.formats.locus_map import render_locus_map_json
from gspec.gbfeatures import FeaturesDocument
from gspec.progress import screen_progress, terminal_rule
from gspec.typing_cluster2 import evaluate_cluster2
from gspec.typing_engine import evaluate_typing
from gspec.typing_models import TypingDocument, uses_cluster_v2, validate_references

_CLUSTER_DEFAULTS = (90.0, 90.0, 96.0)


def reject_gene_engine_flags(
    minid: float, mincov: float, merge_fragments: bool, aligner: AlignerEnum | None
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


def _profile_threshold(cli_value: float, default: float, profile_value: float | None) -> float:
    """Database profile default unless the caller changed the CLI flag.

    The CLI has no explicit "unset" state (its default is the historical
    threshold), so a profile value applies only while the CLI value still
    equals that default. An explicit non-default flag always wins.
    """
    return profile_value if cli_value == default and profile_value is not None else cli_value


def resolve_cluster_params(
    database: Database,
    min_gene_cov: float,
    min_gene_id: float,
    min_cluster_cov: float,
    threads: int,
) -> ClusterParams:
    """Usage-gated constructor: resolve the database profile first, then
    range-check each effective threshold with its CLI flag name before the
    model is built (threads is validated upstream together with the
    gene-path checks)."""
    profile = load_cluster_profile(database.path)
    min_gene_cov = _profile_threshold(min_gene_cov, _CLUSTER_DEFAULTS[0], profile.min_gene_cov)
    min_gene_id = _profile_threshold(min_gene_id, _CLUSTER_DEFAULTS[1], profile.min_gene_id)
    min_cluster_cov = _profile_threshold(
        min_cluster_cov, _CLUSTER_DEFAULTS[2], profile.min_cluster_cov
    )
    for flag, value in (
        ("--min-gene-cov", min_gene_cov),
        ("--min-gene-id", min_gene_id),
        ("--min-cluster-cov", min_cluster_cov),
    ):
        if not 0.0 <= value <= 100.0:
            usage_fail(f"{flag} must be in [0, 100]: got {value}")
    return ClusterParams(
        db=database.name,
        coverage_mode=profile.coverage_mode,
        fallback_evidence=profile.fallback_evidence,
        min_gene_cov=min_gene_cov,
        min_gene_id=min_gene_id,
        min_cluster_cov=min_cluster_cov,
        locus_typeable_coverage=profile.locus_typeable_coverage,
        locus_typeable_identity=profile.locus_typeable_identity,
        locus_quality_min_coverage=(
            profile.locus_quality_min_coverage
            if profile.locus_quality_min_coverage is not None
            else 80.0
        ),
        locus_quality_min_identity=(
            profile.locus_quality_min_identity
            if profile.locus_quality_min_identity is not None
            else 90.0
        ),
        locus_quality_max_missing_genes=(
            profile.locus_quality_max_missing_genes
            if profile.locus_quality_max_missing_genes is not None
            else 3
        ),
        locus_quality_max_extra_panel_genes=(
            profile.locus_quality_max_extra_panel_genes
            if profile.locus_quality_max_extra_panel_genes is not None
            else 2
        ),
        threads=threads,
    )


def _screen_one(
    path: Path,
    database: Database,
    features: FeaturesDocument,
    document: TypingDocument | None,
    params: ClusterParams,
    *,
    quiet: bool,
    debug: bool,
    detail: bool,
) -> ClusterReport:
    """Screen one file and annotate it when the db carries a typing.json —
    one scheme annotates the best call (/1), several schemes populate the
    phenotypes map (/2) — (per-file stderr chatter: Processing / Best
    locus)."""
    if not quiet:
        click.echo(f"Processing: {path}", err=True)
    report = screen_cluster_file(path, database, features, params, debug=debug, detail=detail)
    if document is not None:
        if uses_cluster_v2(document):
            report = report.model_copy(
                update={
                    "phenotypes": evaluate_cluster2(
                        report,
                        document,
                        fallback_evidence=params.fallback_evidence,
                        typeable_coverage=params.locus_typeable_coverage,
                        typeable_identity=params.locus_typeable_identity,
                    )
                }
            )
        else:
            report = evaluate_typing(report, document.schemes[0])
    if report.locus_maps:
        wanted: set[str] = set()
        if report.best is not None:
            wanted.add(report.best.locus)
        wanted.update(
            call.winner_locus
            for call in report.phenotypes.values()
            if call.winner_locus is not None
        )
        report = report.model_copy(
            update={
                "locus_maps": {
                    locus: mapped for locus, mapped in report.locus_maps.items() if locus in wanted
                }
            }
        )
    if not quiet:
        called = report.best.locus if report.best is not None else "none"
        click.echo(f"Best locus in {path}: {called}", err=True)
    return report


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
    detail: bool = False,
    emit: Emit | None = None,
    stdout_emit: bool = False,
    jobs: int = 1,
) -> str:
    """Cluster-engine orchestration beside the gene path: screen each input
    file in order and render (json/md/tsv via formats.cluster*). ``jobs``
    screens N assemblies concurrently, each worker running its own minimap2
    against the shared locus FASTA; stdout still emits in input order
    (head-of-line, file i waits for files 1..i). The db's typing.json is
    reference-validated once up front (TYPING_UNKNOWN_GENE) and picks the
    renderer variant: untyped and single-scheme dbs render gapit.cluster/1
    (byte-frozen), multi-scheme dbs render gapit.cluster/2. With ``emit``,
    chunks stream in document order — tsv/csv/md per completed file (md
    leads with its static frontmatter; each file's section carries its own
    summary row); json is a single document emitted at the end. Chunk
    concatenation is byte-identical to the buffered return value."""
    if jobs < 1:
        usage_fail(f"--jobs must be >= 1: got {jobs}")
    if detail and output_format is not OutputFormat.json:
        usage_fail("--detail requires --format json")
    cpu_count = os.cpu_count()
    if not quiet and cpu_count is not None and jobs * params.threads > cpu_count:
        click.echo(
            f"--jobs {jobs} --threads {params.threads} oversubscribes {cpu_count} cpus", err=True
        )
    features = load_features(database)
    typing_document = load_typing(database)
    v2_document = (
        typing_document
        if typing_document is not None and uses_cluster_v2(typing_document)
        else None
    )
    if typing_document is not None:
        validate_references(typing_document, features)
    typed = typing_document is not None

    def screen_one(path: Path) -> ClusterReport:
        if task_id is not None:
            bar.describe(path.name)
        return _screen_one(
            path,
            database,
            features,
            typing_document,
            params,
            quiet=quiet or task_id is not None,
            debug=debug,
            detail=detail,
        )

    def iter_reports() -> Iterator[ClusterReport]:
        # Each worker owns one minimap2 subprocess. executor.map yields in
        # input order, so rendering and ``emit`` stay head-of-line exactly
        # like the sequential path (a failing file raises at its position).
        def tracked(seq: Iterator[ClusterReport]) -> Iterator[ClusterReport]:
            for report in seq:
                bar.advance()
                yield report

        if jobs == 1:
            yield from tracked(screen_one(path) for path in files)
            return
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            yield from tracked(executor.map(screen_one, files))

    with screen_progress(files, quiet=quiet) as (bar, task_id):
        chunks: list[str] = []

        def sink(chunk: str) -> None:
            if chunk:
                chunks.append(chunk)
                if emit is not None:
                    if stdout_emit:
                        bar.write(chunk)
                    else:
                        with bar.suspend():
                            emit(chunk)

        if output_format is OutputFormat.json:
            reports = list(iter_reports())
            if detail:
                sink(render_locus_map_json(reports, params, now=datetime.now(UTC)))
            elif v2_document is not None:
                sink(render_cluster2_json(reports, params, now=datetime.now(UTC)))
            else:
                sink(render_cluster_json(reports, params, now=datetime.now(UTC), typed=typed))
        elif output_format is OutputFormat.md:
            if v2_document is not None:
                sink(cluster2_md_head(params, now=datetime.now(UTC)))
            else:
                # Static frontmatter leads; each file's section (its own summary
                # row + gene table) streams the moment the file completes.
                sink(cluster_md_head(params, now=datetime.now(UTC)))
            for report in iter_reports():
                if v2_document is not None:
                    sink(cluster2_md_file_chunk(report, v2_document))
                else:
                    sink(cluster_md_file_chunk(report, typed=typed))
        else:
            csv = output_format is OutputFormat.csv
            if v2_document is not None:
                sink(cluster2_tsv_preamble(csv=csv, noheader=noheader, document=v2_document))
                separator = terminal_rule() if stdout_emit and not csv else ""
                if separator and not noheader:
                    sink(separator)
                first = True
                for report in iter_reports():
                    if not first and separator:
                        sink(separator)
                    first = False
                    sink(
                        cluster2_tsv_file_chunk(
                            report, csv=csv, nopath=nopath, document=v2_document
                        )
                    )
                if separator:
                    sink(separator)
            else:
                sink(cluster_tsv_preamble(csv=csv, noheader=noheader, typed=typed))
                separator = terminal_rule() if stdout_emit and not csv else ""
                if separator and not noheader:
                    sink(separator)
                first = True
                for report in iter_reports():
                    if not first and separator:
                        sink(separator)
                    first = False
                    sink(cluster_tsv_file_chunk(report, csv=csv, nopath=nopath, typed=typed))
                if separator:
                    sink(separator)
        return "".join(chunks)
