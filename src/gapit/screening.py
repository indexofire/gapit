"""The blastn contig-screening use-case plus helpers shared by both engines.

The minimap2 engine use-cases (--r1/--r2 reads, --aligner minimap2 assemblies)
live in screening_reads.py; the cluster engine's orchestration and db-kind
guards live in screening_cluster.py (stage-3 split, both under the LOC
ceiling); this module keeps the abricate-parity contig pipeline, the shared
database lookup, and the run_screen dispatcher. OutputFormat/AlignerEnum are
imported from gapit.engines so every screening surface sees the same enums.
"""

import os
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path

import typer

from gapit import db
from gapit.blast import ensure_blast, screen_file
from gapit.bundled import find_bundled, materialize_bundled, resolve_screen_datadir
from gapit.engines import AlignerEnum, Emit, OutputFormat
from gapit.errors import DatabaseError, InputError, ensure_input_file, usage_fail
from gapit.formats.json import render_json
from gapit.formats.md import md_report_chunk, md_report_preamble
from gapit.formats.tsv import tsv_file_chunk, tsv_preamble
from gapit.progress import screen_progress
from gapit.report import Report, ScreeningParams
from gapit.screening_cluster import (
    reject_cluster_engine_flags,
    reject_gene_engine_flags,
    resolve_cluster_params,
    run_cluster_screen,
)


def _resolve_inputs(files: list[Path] | None, fofn: Path | None) -> list[Path]:
    """Input files: --fofn (lines stripped, empties dropped) REPLACES positionals."""
    if fofn is not None:
        if not fofn.is_file():
            raise InputError(
                f"--fofn file not found: {fofn}",
                code="INPUT_NOT_FOUND",
                context={"file": str(fofn)},
            )
        inputs = [
            Path(line.strip())
            for line in fofn.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    elif files:
        inputs = list(files)
    else:
        usage_fail("no input files given (positional FILEs or --fofn)")
    for path in inputs:
        ensure_input_file(path)
    return inputs


def find_database(datadir: Path, name: str, *, quiet: bool = False) -> db.Database:
    """Look up a database by name under the datadir; unknown names list what
    exists.

    A BUNDLED database (gapit.bundled) that is absent from the datadir is
    auto-materialized into it first — deterministic, zero network — so the
    first screen works with no ``db fetch``. Already-materialized bundles
    short-circuit above; a bundle directory left without its certifying
    manifest (interrupted build) is rebuilt by the same path.
    """
    databases = db.discover_databases(datadir)
    for database in databases:
        if database.name == name:
            return database
    if find_bundled(name) is not None and not (datadir / name / "gapit-manifest.json").is_file():
        materialize_bundled(name, datadir, quiet=quiet)
        for database in db.discover_databases(datadir):
            if database.name == name:
                return database
    available = ", ".join(entry.name for entry in databases) or "(none)"
    raise DatabaseError(
        f"Database {name} is not in {datadir}. Available: {available}",
        code="DATABASE_NOT_FOUND",
        context={"db": name, "datadir": str(datadir)},
    )


def run_screen(
    files: list[Path] | None,
    db_name: str,
    datadir: Path | None,
    minid: float,
    mincov: float,
    threads: int,
    jobs: int,
    fofn: Path | None,
    quiet: bool,
    noheader: bool,
    nopath: bool,
    debug: bool,
    output_format: OutputFormat,
    merge_fragments: bool = False,
    aligner: AlignerEnum | None = None,
    min_gene_cov: float = 90.0,
    min_gene_id: float = 90.0,
    min_cluster_cov: float = 96.0,
    emit: Emit | None = None,
) -> str:
    """Screen each input file in order and render; returns the full output
    for the caller to echo (MCP contract), and when ``emit`` is given also
    hands each rendered chunk to it as soon as the chunk exists — tsv/csv/md
    stream per completed file (input order, head-of-line under
    ``--jobs``: file i waits for 1..i; md leads with its static frontmatter,
    counts live in the JSON document only), while json is a single document
    emitted once at the end. Emission order is document order, so chunk
    concatenation is always byte-identical to the buffered return value.

    The per-run gates (blastn presence via ``ensure_blast``, dbtype via one
    ``blastdbcmd -info``) fire once up front, so MISSING_DEPENDENCY and
    DATABASE_NOT_INDEXED errors precede any "Processing:" stderr lines.
    ``merge_fragments`` is the opt-in cross-contig gene-fragment merge
    (gapit extension; default off keeps abricate parity).

    A cluster-kind database dispatches to the minimap2 cluster engine
    (gapit.cluster) after the kind guards: gene-engine flags on a cluster
    db and cluster flags on a gene db are usage errors; the gene path below
    the dispatch is byte-identical to the pre-cluster behavior.
    """
    if not 0.0 < minid <= 100.0:
        usage_fail(f"--minid must be in (0, 100]: got {minid}")
    if not 0.0 <= mincov <= 100.0:
        usage_fail(f"--mincov must be in [0, 100]: got {mincov}")
    if threads < 1:
        usage_fail(f"--threads must be >= 1: got {threads}")
    if jobs < 1:
        usage_fail(f"--jobs must be >= 1: got {jobs}")
    inputs = _resolve_inputs(files, fofn)
    database = find_database(resolve_screen_datadir(datadir, db_name), db_name, quiet=quiet)
    if database.kind == "cluster":
        reject_gene_engine_flags(minid, mincov, merge_fragments, jobs, aligner)
        return run_cluster_screen(
            inputs,
            database,
            resolve_cluster_params(db_name, min_gene_cov, min_gene_id, min_cluster_cov, threads),
            output_format=output_format,
            noheader=noheader,
            nopath=nopath,
            quiet=quiet,
            debug=debug,
            emit=emit,
        )
    reject_cluster_engine_flags(min_gene_cov, min_gene_id, min_cluster_cov)
    params = ScreeningParams(db=db_name, minid=minid, mincov=mincov, threads=threads)
    ensure_blast()
    dbtype = db.blast_db_info(database.sequences_path).dbtype
    cpu_count = os.cpu_count()
    if not quiet and cpu_count is not None and jobs * threads > cpu_count:
        typer.echo(f"--jobs {jobs} --threads {threads} oversubscribes {cpu_count} cpus", err=True)

    def screen_one(path: Path) -> Report:
        if task_id is not None:
            # Bar active: its description carries the current file (notes suppressed).
            bar.describe(path.name)
        elif not quiet:
            typer.echo(f"Processing: {path}", err=True)
        report = screen_file(
            path, database, params, dbtype=dbtype, debug=debug, merge_fragments=merge_fragments
        )
        if task_id is None and not quiet:
            typer.echo(f"Found {len(report.hits)} genes in {path}", err=True)
        return report

    def iter_reports() -> Iterator[Report]:
        # Subprocess-bound work (GIL irrelevant); executor.map yields
        # positionally, so reports stay in input order (SPEC.md §4) and a
        # failing file raises at its position, like the sequential loop.
        def tracked(seq: Iterator[Report]) -> Iterator[Report]:
            for report in seq:  # advance per completed file, input order
                bar.advance()
                yield report

        if jobs == 1:
            yield from tracked(screen_one(path) for path in inputs)
            return
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            yield from tracked(executor.map(screen_one, inputs))

    with screen_progress(inputs, quiet=quiet) as (bar, task_id):
        chunks: list[str] = []

        def sink(chunk: str) -> None:
            if chunk:
                chunks.append(chunk)
                if emit is not None:
                    emit(chunk)

        if output_format is OutputFormat.json:
            sink(render_json(list(iter_reports()), params, now=datetime.now(UTC)))
        elif output_format is OutputFormat.md:
            # Static frontmatter leads (knowable before file 1), then one
            # section per file the moment it completes — exactly like tsv
            # chunks, head-of-line in input order under --jobs.
            sink(md_report_preamble(params, now=datetime.now(UTC)))
            for report in iter_reports():
                sink(md_report_chunk(report))
        else:
            as_csv = output_format is OutputFormat.csv
            sink(tsv_preamble(csv=as_csv, noheader=noheader))
            for report in iter_reports():
                sink(tsv_file_chunk(report, csv=as_csv, nopath=nopath))
    return "".join(chunks)
