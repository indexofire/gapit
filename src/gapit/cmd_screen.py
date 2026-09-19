"""The `gapit screen` command: contig files (blastn) or reads and assemblies
(minimap2) screening.

Lives outside cli.py to keep that module small; cli.py registers it via
``register_screen_command``.
"""

from collections.abc import Callable
from pathlib import Path
from typing import Annotated

import typer

from gapit.errors import GapitError, render_error
from gapit.reads import ReadTypeEnum
from gapit.screening import AlignerEnum, OutputFormat, run_screen, usage_fail
from gapit.screening_reads import run_screen_assemblies, run_screen_reads


def _dispatch(action: Callable[[], None]) -> None:
    """Run a command body; any failure renders the gapit.error/1 envelope on
    stderr and exits with the documented code (mirrors cli._dispatch)."""
    try:
        action()
    except Exception as exc:
        typer.echo(render_error(exc), err=True)
        raise typer.Exit(code=exc.exit_code if isinstance(exc, GapitError) else 1) from exc


Datadir = Annotated[
    Path | None,
    typer.Option(
        "--datadir",
        help="Database directory (default: $GAPIT_DATADIR, then ~/.local/share/gapit/db).",
    ),
]


def screen_command(
    files: Annotated[
        list[Path] | None,
        typer.Argument(help="Input FASTA/GBK/EMBL contig file(s) to screen."),
    ] = None,
    r1: Annotated[
        str | None,
        typer.Option(
            "--r1",
            help="Comma-separated FASTQ reads or assembly FASTA file(s), one per lane.",
        ),
    ] = None,
    r2: Annotated[
        str | None,
        typer.Option("--r2", help="Comma-separated mate FASTQ file(s); must match --r1 count."),
    ] = None,
    read_type: Annotated[
        ReadTypeEnum | None,
        typer.Option(
            "--read-type",
            help=(
                "minimap2 preset for reads mode (default: sr for FASTQ, map-ont"
                " for assembly FASTA)."
            ),
        ),
    ] = None,
    min_breadth: Annotated[
        float,
        typer.Option("--min-breadth", help="Reads mode: minimum %breadth for presence."),
    ] = 90.0,
    aligner: Annotated[
        AlignerEnum | None,
        typer.Option(
            "--aligner",
            help=(
                "Alignment engine (default: blastn for contig files, minimap2 for --r1/--r2 reads)."
            ),
        ),
    ] = None,
    db: Annotated[
        str, typer.Option("--db", help="Database to screen against (datadir subdir).")
    ] = "ncbi",
    datadir: Datadir = None,
    minid: Annotated[
        float, typer.Option("--minid", help="Minimum %identity, 0 < x <= 100.")
    ] = 80.0,
    mincov: Annotated[
        float, typer.Option("--mincov", help="Minimum %coverage, 0 <= x <= 100.")
    ] = 80.0,
    threads: Annotated[int, typer.Option("--threads", help="BLAST worker threads.")] = 1,
    jobs: Annotated[
        int,
        typer.Option(
            "--jobs",
            help=(
                "Screen N input files concurrently (gapit extension; output order is"
                " always input order). Each worker runs its own BLAST against the"
                " shared db index, which BLAST mmaps — concurrent readers are fine."
            ),
        ),
    ] = 1,
    fofn: Annotated[
        Path | None,
        typer.Option("--fofn", help="File of filenames; replaces the positional FILEs."),
    ] = None,
    quiet: Annotated[bool, typer.Option("--quiet", help="Silence stderr diagnostics.")] = False,
    noheader: Annotated[bool, typer.Option("--noheader", help="Suppress the header row.")] = False,
    nopath: Annotated[bool, typer.Option("--nopath", help="Basename the FILE column.")] = False,
    debug: Annotated[bool, typer.Option("--debug", help="Verbose stderr diagnostics.")] = False,
    output_format: Annotated[
        OutputFormat | None,
        typer.Option("--format", help="Output format (reads mode defaults to json)."),
    ] = None,
) -> None:
    """Screen contig files or FASTQ reads (R1 and R2 comma-lists, one lane
    each) for AMR/virulence genes."""

    def run() -> None:
        if (r1 is not None or r2 is not None) and files:
            usage_fail("--r1/--r2 and positional contig FILEs are mutually exclusive")
        if r2 is not None and r1 is None:
            usage_fail("--r2 requires --r1")
        if r1 is not None or r2 is not None:
            run_screen_reads(
                r1 or "",
                r2,
                db,
                datadir,
                read_type,
                min_breadth,
                threads,
                output_format,
                quiet,
                debug,
                aligner=aligner,
            )
        elif aligner is AlignerEnum.minimap2:
            run_screen_assemblies(
                files,
                fofn,
                db,
                datadir,
                read_type,
                min_breadth,
                threads,
                output_format,
                quiet,
                debug,
            )
        else:
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
            )

    _dispatch(run)


def register_screen_command(app: typer.Typer) -> None:
    """Attach the screen command to the CLI app."""
    app.command("screen")(screen_command)
