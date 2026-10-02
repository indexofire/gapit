"""The `gapit screen` command: contig files (blastn) or reads and assemblies
(minimap2) screening.

Lives outside cli.py to keep that module small; cli.py registers it via
``register_screen_command``. The reads-mode argument rules (comma-list
splitting, positional-FASTQ classification, flag rejects) live in
cmd_screen_reads_args.py.
"""

from pathlib import Path
from typing import Annotated

import typer

from gapit.cmd_screen_positionals import reject_db_name_positionals
from gapit.cmd_screen_reads_args import (
    reads_positional,
    reject_mixed_positionals,
    reject_reads_mode_flags,
    split_read_list,
)
from gapit.dispatch import Datadir, dispatch, output_target
from gapit.errors import usage_fail
from gapit.reads import ReadTypeEnum
from gapit.screening import (
    AlignerEnum,
    OutputFormat,
    run_screen,
)
from gapit.screening_cluster import reject_cluster_engine_flags
from gapit.screening_reads import run_screen_assemblies, run_screen_reads
from gapit.screening_reads_positional import run_screen_reads_positional


def screen_command(
    files: Annotated[
        list[Path] | None,
        typer.Argument(
            help=(
                "Input contig file(s) to screen; an all-FASTQ wildcard enters reads"
                " mode with samples auto-paired from filenames."
            ),
        ),
    ] = None,
    r1: Annotated[
        str | None,
        typer.Option(
            "--r1",
            "-1",
            help="Comma-separated FASTQ reads or assembly FASTA file(s), one per lane.",
        ),
    ] = None,
    r2: Annotated[
        str | None,
        typer.Option(
            "--r2", "-2", help="Comma-separated mate FASTQ file(s); must match --r1 count."
        ),
    ] = None,
    read_type: Annotated[
        ReadTypeEnum | None,
        typer.Option(
            "--read-type",
            "-x",
            help=(
                "minimap2 preset for reads mode (default: sr for FASTQ, map-ont"
                " for assembly FASTA)."
            ),
        ),
    ] = None,
    min_breadth: Annotated[
        float,
        typer.Option("--min-breadth", "-b", help="Reads mode: minimum %breadth for presence."),
    ] = 90.0,
    min_identity: Annotated[
        float,
        typer.Option(
            "--min-identity",
            "-I",
            help=(
                "Reads mode: minimum %identity per alignment, 0 <= x <= 100 (0 = off;"
                " any nonzero value emits gapit.reads/2)."
            ),
        ),
    ] = 0.0,
    min_mapq: Annotated[
        int,
        typer.Option(
            "--min-mapq",
            "-M",
            help=(
                "Reads mode: minimum MAPQ per alignment (0 = off; any nonzero value"
                " emits gapit.reads/2)."
            ),
        ),
    ] = 0,
    aligner: Annotated[
        AlignerEnum | None,
        typer.Option(
            "--aligner",
            "-a",
            help=(
                "Alignment engine (default: blastn for contig files, minimap2 for --r1/--r2 reads)."
            ),
        ),
    ] = None,
    datadir: Datadir = None,
    minid: Annotated[
        float, typer.Option("--minid", "-i", help="Minimum %identity, 0 < x <= 100.")
    ] = 80.0,
    mincov: Annotated[
        float, typer.Option("--mincov", "-c", help="Minimum %coverage, 0 <= x <= 100.")
    ] = 80.0,
    min_gene_cov: Annotated[
        float,
        typer.Option(
            "--min-gene-cov", "-g", help="Cluster dbs: min %coverage for a present gene verdict."
        ),
    ] = 90.0,
    min_gene_id: Annotated[
        float,
        typer.Option(
            "--min-gene-id", "-G", help="Cluster dbs: min %identity for a present gene verdict."
        ),
    ] = 90.0,
    min_cluster_cov: Annotated[
        float,
        typer.Option(
            "--min-cluster-cov", "-C", help="Cluster dbs: min locus %coverage for a best call."
        ),
    ] = 96.0,
    threads: Annotated[int, typer.Option("--threads", "-t", help="BLAST worker threads.")] = 1,
    jobs: Annotated[
        int,
        typer.Option(
            "--jobs",
            "-j",
            help=(
                "Screen N input files concurrently (contig files, or the wildcard"
                " FASTQ path's samples; gapit extension; output order is always"
                " input order). Each worker runs its own aligner against the"
                " shared db — BLAST mmaps the index and minimap2 holds it in"
                " memory, so concurrent readers are fine."
            ),
        ),
    ] = 1,
    all_genes: Annotated[
        bool,
        typer.Option(
            "--all-genes",
            "-A",
            help="Reads table: also list absent gene calls (default: present genes only).",
        ),
    ] = False,
    merge_fragments: Annotated[
        bool,
        typer.Option(
            "--merge-fragments/--no-merge-fragments",
            "-m",
            help=(
                "Merge gene fragments split across contigs (gapit extension): report one"
                " hit when fragments of a gene that each fail --mincov jointly cover >= mincov"
                " of the subject. blastn contig mode only."
            ),
        ),
    ] = False,
    fofn: Annotated[
        Path | None,
        typer.Option("--fofn", "-F", help="File of filenames; replaces the positional FILEs."),
    ] = None,
    quiet: Annotated[
        bool, typer.Option("--quiet", "-q", help="Silence stderr diagnostics.")
    ] = False,
    noheader: Annotated[
        bool, typer.Option("--noheader", "-n", help="Suppress the header row.")
    ] = False,
    nopath: Annotated[
        bool, typer.Option("--nopath", "-p", help="Basename the FILE column.")
    ] = False,
    debug: Annotated[
        bool, typer.Option("--debug", "-v", help="Verbose stderr diagnostics.")
    ] = False,
    output_format: Annotated[
        OutputFormat | None,
        typer.Option(
            "--format",
            "-f",
            help=(
                "Output format (tsv is the default everywhere: reads mode"
                " streams the table per completed file or sample). json/md are"
                " the explicit agent opt-ins; json is written once at the end"
                " (single document)."
            ),
        ),
    ] = None,
    output: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help=(
                "Write the report to PATH instead of stdout (truncates; streamed chunks"
                " are flushed per file; stdout then carries no data)."
            ),
        ),
    ] = None,
    *,
    db: Annotated[
        str,
        typer.Option(
            "--db",
            "-d",
            "-db",
            help="Database to screen against (required, no default).",
        ),
    ],
) -> None:
    """Screen contig files or FASTQ reads (R1 and R2 comma-lists, one lane
    each, or an all-FASTQ positional wildcard auto-paired into samples) for
    known genes (reference or custom databases)."""

    def run() -> None:
        args = list(files or [])
        if (r1 is not None or r2 is not None) and args:
            usage_fail("--r1/--r2 and positional contig FILEs are mutually exclusive")
        if r2 is not None and r1 is None:
            usage_fail("--r2 requires --r1")
        # Auto reads mode only when the default engine would otherwise run:
        # an explicit --aligner keeps its own routing (minimap2 assemblies,
        # blastn contigs) and the frozen matrix of guards with it.
        reads_flags = [reads_positional(path) for path in args] if args and aligner is None else []
        all_reads = bool(reads_flags) and all(reads_flags)
        if (
            (min_identity > 0 or min_mapq > 0)
            and r1 is None
            and aligner is not AlignerEnum.minimap2
            and not all_reads
        ):
            usage_fail("--min-identity/--min-mapq are reads-mode only (minimap2 engine)")
        reject_db_name_positionals(files, datadir)
        with output_target(output) as deliver:
            if r1 is not None or r2 is not None:
                reject_reads_mode_flags(fofn, noheader, nopath, jobs, merge_fragments)
                reject_cluster_engine_flags(min_gene_cov, min_gene_id, min_cluster_cov)
                deliver(
                    run_screen_reads(
                        split_read_list(r1 or "", "--r1"),
                        split_read_list(r2, "--r2") if r2 is not None else None,
                        db,
                        datadir,
                        read_type,
                        min_breadth,
                        min_identity,
                        min_mapq,
                        threads,
                        output_format,
                        quiet,
                        debug,
                        aligner=aligner,
                        minid=minid,
                        mincov=mincov,
                        all_genes=all_genes,
                    )
                )
            elif any(reads_flags):
                reject_mixed_positionals(args, reads_flags)
                # --jobs is legal ONLY here in reads mode (the wildcard
                # screens many samples); passing 1 to the shared guard keeps
                # the other reads-mode rejects frozen, and the use-case owns
                # the jobs >= 1 validation like run_screen.
                reject_reads_mode_flags(fofn, noheader, nopath, 1, merge_fragments)
                reject_cluster_engine_flags(min_gene_cov, min_gene_id, min_cluster_cov)

                def warn(message: str) -> None:
                    if not quiet:
                        typer.echo(f"WARNING: {message}", err=True)

                # run_screen_reads_positional emits its output through
                # `deliver` (tsv/md per completed sample, json once); echoing
                # the return value here would duplicate every byte.
                run_screen_reads_positional(
                    args,
                    db,
                    datadir,
                    read_type,
                    min_breadth,
                    min_identity,
                    min_mapq,
                    threads,
                    output_format,
                    quiet,
                    warn,
                    debug,
                    minid=minid,
                    mincov=mincov,
                    jobs=jobs,
                    emit=deliver,
                    all_genes=all_genes,
                )
            elif aligner is AlignerEnum.minimap2:
                if merge_fragments:
                    usage_fail("--merge-fragments is not available with --aligner minimap2")
                reject_cluster_engine_flags(min_gene_cov, min_gene_id, min_cluster_cov)
                deliver(
                    run_screen_assemblies(
                        files,
                        fofn,
                        db,
                        datadir,
                        read_type,
                        min_breadth,
                        min_identity,
                        min_mapq,
                        threads,
                        jobs,
                        noheader,
                        nopath,
                        output_format,
                        quiet,
                        debug,
                        minid=minid,
                        mincov=mincov,
                    )
                )
            else:
                # run_screen emits its output through `deliver`; echoing the
                # return value here would duplicate every byte.
                run_screen(
                    files,
                    db,
                    datadir,
                    minid,
                    mincov,
                    threads,
                    jobs,
                    fofn,
                    quiet,
                    noheader,
                    nopath,
                    debug,
                    output_format or OutputFormat.tsv,
                    merge_fragments=merge_fragments,
                    aligner=aligner,
                    min_gene_cov=min_gene_cov,
                    min_gene_id=min_gene_id,
                    min_cluster_cov=min_cluster_cov,
                    emit=deliver,
                )

    dispatch(run)


def register_screen_command(app: typer.Typer) -> None:
    """Attach the screen command to the CLI app (bare invocation prints help)."""
    app.command("screen", no_args_is_help=True)(screen_command)
