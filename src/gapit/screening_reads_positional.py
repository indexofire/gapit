"""The positional-FASTQ wildcard screening use-case (SPEC.md §10 extension).

``gapit screen -d db *.fastq.gz``: when every positional file is FASTQ (by
extension, or content sniff for ambiguous names — the classification lives
in cmd_screen_reads_args.reads_positional), the command enters reads mode
without ``--r1``/``--r2`` and gapit.readpairs infers samples from filename
conventions. The document carries ONE files[] entry per sample (the union
over that sample's lanes) keyed by the sample name; the single-sample
``--r1``/``--r2`` path in screening_reads.py is untouched.

Streaming/parallel contract (the run_screen mirror): ``--jobs N`` screens
samples concurrently — executor.map yields positionally, so results stay in
sample order and a failing sample raises at its position, like the
sequential loop. With ``emit`` the streaming formats lead with a static
preamble first and one chunk per completed sample, head-of-line in sample
order under ``--jobs``; json stays a single document emitted once at the
end. Chunk concatenation is always byte-identical to the buffered return
value.
"""

import os
from collections.abc import Callable, Iterator
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path

import typer

from gapit.bundled import resolve_screen_datadir
from gapit.engines import Emit
from gapit.errors import ensure_input_file, usage_fail
from gapit.formats.md import reads_md_chunk, reads_md_preamble
from gapit.formats.reads_json import render_reads2_json, render_reads_json
from gapit.formats.reads_tsv import reads_tsv_chunk, reads_tsv_preamble
from gapit.progress import screen_progress
from gapit.readpairs import SampleLanes, pair_samples
from gapit.reads import ReadsParams, ReadsReport, ReadTypeEnum, screen_reads
from gapit.screening import OutputFormat, find_database
from gapit.screening_reads import (
    reject_blastn_thresholds,
    resolve_read_preset,
    validate_reads_usage,
)


def run_screen_reads_positional(
    files: list[Path],
    db_name: str,
    datadir: Path | None,
    read_type: ReadTypeEnum | None,
    min_breadth: float,
    min_identity: float,
    min_mapq: int,
    threads: int,
    output_format: OutputFormat | None,
    quiet: bool,
    warn: Callable[[str], None],
    debug: bool = False,
    minid: float = 80.0,
    mincov: float = 80.0,
    jobs: int = 1,
    emit: Emit | None = None,
    all_genes: bool = False,
) -> str:
    """Screen positional FASTQ file(s) as auto-paired samples: one minimap2
    run per lane, per-sample union, one files[] entry per sample (reads[]
    carries the sample key; the lane files parade on stderr). tsv is the
    default format (streaming: header, then one chunk per sample); json
    stays the opt-in single document (SPEC.md §10); a nonzero
    --min-identity/--min-mapq turns on gapit.reads/2 alignment filtering.
    Returns the full output for the caller to echo (MCP contract), and when
    ``emit`` is given also hands each rendered chunk to it as soon as the
    chunk exists (module docstring)."""
    validate_reads_usage(min_breadth, min_identity, min_mapq, threads)
    reject_blastn_thresholds(minid, mincov, "positional reads")
    if jobs < 1:
        usage_fail(f"--jobs must be >= 1: got {jobs}")
    for path in files:
        ensure_input_file(path, "reads file")
    database = find_database(resolve_screen_datadir(datadir, db_name), db_name, quiet=quiet)
    if database.kind == "cluster":
        usage_fail(
            "cluster databases are assembly-contig screening only;"
            " positional FASTQ reads are not available"
        )
    samples = pair_samples(files, warn=warn)
    lanes = [lane for sample in samples for lane in sample.lanes]
    resolved = resolve_read_preset(lanes, read_type, quiet)
    cpu_count = os.cpu_count()
    if not quiet and cpu_count is not None and jobs * threads > cpu_count:
        typer.echo(f"--jobs {jobs} --threads {threads} oversubscribes {cpu_count} cpus", err=True)

    def screen_sample(sample: SampleLanes) -> ReadsReport:
        if task_id is not None:
            # Bar active: its description carries the sample (notes suppressed).
            bar.describe(sample.sample)
        elif not quiet:
            read_list = ", ".join(
                str(path) for lane in sample.lanes for path in lane if path is not None
            )
            typer.echo(f"Screening sample {sample.sample} reads: {read_list}", err=True)
        report = screen_reads(
            list(sample.lanes),
            database,
            read_type=resolved.value,
            min_breadth=min_breadth,
            threads=threads,
            debug=debug,
            min_identity=min_identity,
            min_mapq=min_mapq,
        )
        present = sum(1 for gene in report.genes if gene.present)
        if task_id is None and not quiet:
            typer.echo(f"Detected {present} present genes in sample {sample.sample}", err=True)
        return report.model_copy(update={"reads": (sample.sample,)})

    def iter_reports() -> Iterator[ReadsReport]:
        # Subprocess-bound work (GIL irrelevant); executor.map yields
        # positionally, so reports stay in sample order and a failing sample
        # raises at its position, like the sequential loop (run_screen twin).
        def tracked(seq: Iterator[ReadsReport]) -> Iterator[ReadsReport]:
            for report in seq:  # advance per completed sample, sample order
                bar.advance()
                yield report

        if jobs == 1:
            yield from tracked(screen_sample(sample) for sample in samples)
            return
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            yield from tracked(executor.map(screen_sample, samples))

    params = ReadsParams(
        db=db_name,
        read_type=resolved.value,
        min_breadth=min_breadth,
        threads=threads,
        min_identity=min_identity,
        min_mapq=min_mapq,
    )
    now = datetime.now(UTC)
    reads2 = min_identity > 0.0 or min_mapq > 0
    fmt = output_format if output_format is not None else OutputFormat.tsv
    with screen_progress(samples, quiet=quiet) as (bar, task_id):
        chunks: list[str] = []

        def sink(chunk: str) -> None:
            if chunk:
                chunks.append(chunk)
                if emit is not None:
                    emit(chunk)

        if fmt is OutputFormat.md:
            # Static frontmatter leads (knowable before sample 1), then one
            # section per sample the moment it completes — head-of-line in
            # sample order under --jobs (the cluster_md precedent; run totals
            # live in the JSON document only).
            sink(reads_md_preamble(params, now=now, reads2=reads2))
            for report in iter_reports():
                sink(reads_md_chunk(report, reads2=reads2))
        elif fmt is OutputFormat.tsv or fmt is OutputFormat.csv:
            # The default table format: the header is emitted LAZILY with the
            # first sample's chunk — per-sample stderr progress ("Screening
            # sample …", "Detected … genes") then precedes the table on the
            # terminal (rightsholder UX), instead of the header printing
            # before any sample has started.
            csv = fmt is OutputFormat.csv
            header_pending = True
            for report in iter_reports():
                if header_pending:
                    sink(reads_tsv_preamble(reads2=reads2, csv=csv))
                    header_pending = False
                sink(reads_tsv_chunk(report, reads2=reads2, csv=csv, all_genes=all_genes))
            if header_pending and not chunks:
                sink(reads_tsv_preamble(reads2=reads2, csv=csv))
        else:
            render = render_reads2_json if reads2 else render_reads_json
            sink(render(list(iter_reports()), params, now=now))
    return "".join(chunks)
