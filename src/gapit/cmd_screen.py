"""The `gapit screen` command: contig files (blastn) or reads and assemblies
(minimap2) screening.

Lives outside cli.py to keep that module small; cli.py registers it via
``register_screen_command``. The reads-mode argument rules (comma-list
splitting, positional-FASTQ classification, flag rejects) live in
cmd_screen_reads_args.py.
"""

from pathlib import Path
from typing import TypeVar

import click
from click.core import ParameterSource

from gapit.cmd_screen_positionals import reject_db_name_positionals
from gapit.cmd_screen_reads_args import (
    reads_positional,
    reject_mixed_positionals,
    reject_reads_mode_flags,
    split_read_list,
)
from gapit.dispatch import EnumChoice, datadir_option, dispatch, output_target, quiet_option
from gapit.errors import usage_fail
from gapit.reads import ReadTypeEnum
from gapit.screen_defaults import screen_defaults_for
from gapit.screening import (
    AlignerEnum,
    OutputFormat,
    run_screen,
)
from gapit.screening_cluster import reject_cluster_engine_flags
from gapit.screening_reads import run_screen_assemblies, run_screen_reads
from gapit.screening_reads_positional import run_screen_reads_positional


def _reject_reads_only_flags(
    read_type: ReadTypeEnum | None,
    min_breadth: float,
    all_genes: bool,
    min_depth: float = 0.0,
    min_relative_depth: float = 0.0,
    relative_to: str = "lpxA",
) -> None:
    """Reads-engine flags are usage errors outside reads mode, never silent."""
    if read_type is not None:
        usage_fail(
            "--read-type is reads-mode only (--r1/--r2, positional FASTQ, or --aligner minimap2)"
        )
    if min_breadth != 90.0:
        usage_fail(
            "--min-breadth is reads-mode only (--r1/--r2, positional FASTQ, or --aligner minimap2)"
        )
    if min_depth != 0.0:
        usage_fail(
            "--min-depth is reads-mode only (--r1/--r2, positional FASTQ, or --aligner minimap2)"
        )
    if min_relative_depth != 0.0:
        usage_fail(
            "--min-relative-depth is reads-mode only"
            " (--r1/--r2, positional FASTQ, or --aligner minimap2)"
        )
    if relative_to != "lpxA":
        usage_fail(
            "--relative-to is reads-mode only (--r1/--r2, positional FASTQ, or --aligner minimap2)"
        )
    if all_genes:
        usage_fail(
            "--all-genes is reads-mode only (--r1/--r2, positional FASTQ, or --aligner minimap2)"
        )


@click.argument(
    "files",
    nargs=-1,
    required=False,
    type=click.Path(path_type=Path),
    help=(
        "Input contig file(s) to screen; an all-FASTQ wildcard enters reads"
        " mode with samples auto-paired from filenames."
    ),
)
@click.option(
    "--r1",
    "-1",
    help="Comma-separated FASTQ reads or assembly FASTA file(s), one per lane.",
)
@click.option("--r2", "-2", help="Comma-separated mate FASTQ file(s); must match --r1 count.")
@click.option(
    "--read-type",
    "-x",
    type=EnumChoice(ReadTypeEnum),
    help=("minimap2 preset for reads mode (default: sr for FASTQ, map-ont for assembly FASTA)."),
)
@click.option(
    "--min-breadth",
    "-b",
    default=90.0,
    help="Reads mode: minimum %breadth for presence.",
)
@click.option(
    "--min-identity",
    "-I",
    default=0.0,
    help=(
        "Reads mode: minimum %identity per alignment, 0 <= x <= 100 (0 = off;"
        " any nonzero value emits gapit.reads/2)."
    ),
)
@click.option(
    "--min-mapq",
    "-M",
    default=0,
    help=(
        "Reads mode: minimum MAPQ per alignment (0 = off; any nonzero value emits gapit.reads/2)."
    ),
)
@click.option(
    "--min-depth",
    "-e",
    default=0.0,
    help="Reads mode: minimum mean depth for presence (0 = off).",
)
@click.option(
    "--min-relative-depth",
    "-R",
    default=0.0,
    help=(
        "Reads mode: minimum gene:control mean-depth ratio (0..1; 0 = off)."
        " The control gene is --relative-to."
    ),
)
@click.option(
    "--relative-to",
    "-r",
    default="lpxA",
    help="Reads mode: control gene for --min-relative-depth (default: lpxA).",
)
@click.option(
    "--aligner",
    "-a",
    type=EnumChoice(AlignerEnum),
    help=("Alignment engine (default: blastn for contig files, minimap2 for --r1/--r2 reads)."),
)
@datadir_option()
@click.option("--minid", "-i", default=80.0, help="Minimum %identity, 0 < x <= 100.")
@click.option("--mincov", "-c", default=80.0, help="Minimum %coverage, 0 <= x <= 100.")
@click.option(
    "--min-gene-cov",
    "-g",
    default=90.0,
    help="Cluster dbs: min %coverage for a present gene verdict.",
)
@click.option(
    "--min-gene-id",
    "-G",
    default=90.0,
    help="Cluster dbs: min %identity for a present gene verdict.",
)
@click.option(
    "--min-cluster-cov",
    "-C",
    default=96.0,
    help="Cluster dbs: min locus %coverage for a best call.",
)
@click.option("--threads", "-t", default=1, help="BLAST worker threads.")
@click.option(
    "--jobs",
    "-j",
    default=1,
    help=(
        "Screen N input files concurrently (contig files, cluster-db"
        " assemblies, or the wildcard FASTQ path's samples; gapit"
        " extension; output order is always input order). Each worker runs"
        " its own aligner against the shared db — BLAST mmaps the index and"
        " minimap2 holds it in memory, so concurrent readers are fine."
    ),
)
@click.option(
    "--all-genes",
    "-A",
    is_flag=True,
    help="Reads table: also list absent gene calls (default: present genes only).",
)
@click.option(
    "--merge-fragments/--no-merge-fragments",
    "-m",
    is_flag=True,
    help=(
        "Merge gene fragments split across contigs (gapit extension): report one"
        " hit when fragments of a gene that each fail --mincov jointly cover >= mincov"
        " of the subject. blastn contig mode only."
    ),
)
@click.option(
    "--fofn",
    "-F",
    type=click.Path(path_type=Path),
    help="File of filenames; replaces the positional FILEs.",
)
@quiet_option()
@click.option("--noheader", "-n", is_flag=True, help="Suppress the header row.")
@click.option("--nopath", "-p", is_flag=True, help="Basename the FILE column.")
@click.option("--debug", "-v", is_flag=True, help="Verbose stderr diagnostics.")
@click.option(
    "--format",
    "-f",
    "output_format",
    type=EnumChoice(OutputFormat),
    help=(
        "Output format (tsv is the default everywhere: reads mode"
        " streams the table per completed file or sample). json/md are"
        " the explicit agent opt-ins; json is written once at the end"
        " (single document)."
    ),
)
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    help=(
        "Write the report to PATH instead of stdout (truncates; streamed chunks"
        " are flushed per file; stdout then carries no data)."
    ),
)
@click.option(
    "--db",
    "-d",
    "-db",
    required=True,
    help="Database to screen against (required, no default).",
)
def screen_command(
    files: tuple[Path, ...],
    r1: str | None,
    r2: str | None,
    read_type: ReadTypeEnum | None,
    min_breadth: float,
    min_identity: float,
    min_mapq: int,
    min_depth: float,
    min_relative_depth: float,
    relative_to: str,
    aligner: AlignerEnum | None,
    datadir: Path | None,
    minid: float,
    mincov: float,
    min_gene_cov: float,
    min_gene_id: float,
    min_cluster_cov: float,
    threads: int,
    jobs: int,
    all_genes: bool,
    merge_fragments: bool,
    fofn: Path | None,
    quiet: bool,
    noheader: bool,
    nopath: bool,
    debug: bool,
    output_format: OutputFormat | None,
    output: Path | None,
    *,
    db: str,
) -> None:
    """Screen contig files or FASTQ reads (R1 and R2 comma-lists, one lane
    each, or an all-FASTQ positional wildcard auto-paired into samples) for
    known genes (reference or custom databases)."""
    file_args: list[Path] | None = list(files) if files else None

    def run() -> None:
        args = list(file_args or [])
        if (r1 is not None or r2 is not None) and args:
            usage_fail("--r1/--r2 and positional contig FILEs are mutually exclusive")
        if r2 is not None and r1 is None:
            usage_fail("--r2 requires --r1")
        # Auto reads mode only when the default engine would otherwise run:
        # an explicit --aligner keeps its own routing (minimap2 assemblies,
        # blastn contigs) and the frozen matrix of guards with it.
        reads_flags = [reads_positional(path) for path in args] if args and aligner is None else []
        all_reads = bool(reads_flags) and all(reads_flags)
        defaults = screen_defaults_for(db, datadir)
        ctx = click.get_current_context()

        T = TypeVar("T")

        def resolve(param: str, current: T, db_value: T | None) -> T:
            """CLI value wins; otherwise the database default, then the CLI
            default. Click records DEFAULT for options the user did not pass."""
            if ctx.get_parameter_source(param) is ParameterSource.COMMANDLINE:
                return current
            return current if db_value is None else db_value

        reads_db = defaults.reads if defaults is not None else None
        contigs_db = defaults.contigs if defaults is not None else None
        eff_breadth = resolve(
            "min_breadth", min_breadth, reads_db.min_breadth if reads_db else None
        )
        eff_identity = resolve(
            "min_identity", min_identity, reads_db.min_identity if reads_db else None
        )
        eff_mapq = resolve("min_mapq", min_mapq, reads_db.min_mapq if reads_db else None)
        eff_depth = resolve("min_depth", min_depth, reads_db.min_depth if reads_db else None)
        eff_relative = resolve(
            "min_relative_depth",
            min_relative_depth,
            reads_db.min_relative_depth if reads_db else None,
        )
        eff_relative_to = resolve(
            "relative_to", relative_to, reads_db.relative_to if reads_db else None
        )
        eff_minid = resolve("minid", minid, contigs_db.minid if contigs_db else None)
        eff_mincov = resolve("mincov", mincov, contigs_db.mincov if contigs_db else None)
        eff_merge = resolve(
            "merge_fragments",
            merge_fragments,
            contigs_db.merge_fragments if contigs_db else None,
        )
        if (
            (min_identity > 0 or min_mapq > 0)
            and r1 is None
            and aligner is not AlignerEnum.minimap2
            and not all_reads
        ):
            usage_fail("--min-identity/--min-mapq are reads-mode only (minimap2 engine)")
        reject_db_name_positionals(file_args, datadir)
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
                        eff_breadth,
                        eff_identity,
                        eff_mapq,
                        threads,
                        output_format,
                        quiet,
                        debug,
                        aligner=aligner,
                        minid=minid,
                        mincov=mincov,
                        all_genes=all_genes,
                        min_depth=eff_depth,
                        min_relative_depth=eff_relative,
                        relative_to=eff_relative_to,
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
                        click.echo(f"WARNING: {message}", err=True)

                # run_screen_reads_positional emits its output through
                # `deliver` (tsv/md per completed sample, json once); echoing
                # the return value here would duplicate every byte.
                run_screen_reads_positional(
                    args,
                    db,
                    datadir,
                    read_type,
                    eff_breadth,
                    eff_identity,
                    eff_mapq,
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
                    stdout_emit=output is None,
                    min_depth=eff_depth,
                    min_relative_depth=eff_relative,
                    relative_to=eff_relative_to,
                )
            elif aligner is AlignerEnum.minimap2:
                if merge_fragments:
                    usage_fail("--merge-fragments is not available with --aligner minimap2")
                reject_cluster_engine_flags(min_gene_cov, min_gene_id, min_cluster_cov)
                deliver(
                    run_screen_assemblies(
                        file_args,
                        fofn,
                        db,
                        datadir,
                        read_type,
                        eff_breadth,
                        eff_identity,
                        eff_mapq,
                        threads,
                        jobs,
                        noheader,
                        nopath,
                        output_format,
                        quiet,
                        debug,
                        minid=minid,
                        mincov=mincov,
                        all_genes=all_genes,
                        min_depth=eff_depth,
                        min_relative_depth=eff_relative,
                        relative_to=eff_relative_to,
                    )
                )
            else:
                _reject_reads_only_flags(
                    read_type,
                    min_breadth,
                    all_genes,
                    min_depth,
                    min_relative_depth,
                    relative_to,
                )
                # run_screen emits its output through `deliver`; echoing
                # the return value here would duplicate every byte.
                run_screen(
                    file_args,
                    db,
                    datadir,
                    eff_minid,
                    eff_mincov,
                    threads,
                    jobs,
                    fofn,
                    quiet,
                    noheader,
                    nopath,
                    debug,
                    output_format or OutputFormat.tsv,
                    merge_fragments=eff_merge,
                    aligner=aligner,
                    min_gene_cov=min_gene_cov,
                    min_gene_id=min_gene_id,
                    min_cluster_cov=min_cluster_cov,
                    emit=deliver,
                    stdout_emit=output is None,
                )

    dispatch(run)


def register_screen_command(app: click.Group) -> None:
    """Attach the screen command to the CLI app (bare invocation prints help)."""
    app.command("screen", no_args_is_help=True)(screen_command)
