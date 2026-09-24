"""The minimap2 screening use-cases (SPEC.md §10): ``--r1``/``--r2`` reads and
positional assemblies via ``--aligner minimap2``.

Split from screening.py so each use-case module stays under the 250 pure-LOC
ceiling; screening.py keeps the blastn contig pipeline and shared helpers.
"""

from datetime import UTC, datetime
from pathlib import Path

import typer

from gapit import config
from gapit.errors import ensure_input_file, usage_fail
from gapit.formats.json import render_reads2_json, render_reads_json
from gapit.formats.md import render_reads2_markdown, render_reads_markdown
from gapit.reads import (
    ReadFileKind,
    ReadsParams,
    ReadTypeEnum,
    detect_read_kind,
    screen_reads,
)
from gapit.screening import AlignerEnum, OutputFormat, find_database


def _pair_read_lanes(r1: list[Path], r2: list[Path] | None) -> list[tuple[Path, Path | None]]:
    """Pair per-lane --r1/--r2 file lists (r2 None = single-end)."""
    if r2 is not None and len(r2) != len(r1):
        usage_fail(f"--r2 has {len(r2)} files but --r1 has {len(r1)} (lanes must pair up)")
    if r2 is None:
        return [(path, None) for path in r1]
    return list(zip(r1, r2, strict=True))


def _resolve_read_preset(
    lanes: list[tuple[Path, Path | None]], read_type: ReadTypeEnum | None, quiet: bool
) -> ReadTypeEnum:
    """Detect every input's kind from content and resolve the minimap2 preset:
    assembly FASTA forces map-ont (one stderr note unless quiet/explicit),
    FASTQ keeps sr unless a preset was given. Mixed kinds and FASTA paired-end
    are usage errors."""
    r1_kinds = {detect_read_kind(r1_path) for r1_path, _ in lanes}
    if len(r1_kinds) > 1:
        usage_fail("mixed FASTA and FASTQ inputs")
    r2_kinds = {detect_read_kind(r2_path) for _, r2_path in lanes if r2_path is not None}
    if r2_kinds and ReadFileKind.fasta in r1_kinds | r2_kinds:
        usage_fail("paired-end requires FASTQ")
    if ReadFileKind.fasta not in r1_kinds:
        return read_type if read_type is not None else ReadTypeEnum.sr
    if read_type is not None and read_type is not ReadTypeEnum.map_ont:
        usage_fail("assembly FASTA requires map-ont")
    if read_type is None and not quiet:
        typer.echo("assembly FASTA detected; using map-ont", err=True)
    return ReadTypeEnum.map_ont


def _reject_blastn_thresholds(minid: float, mincov: float, mode: str) -> None:
    """--minid/--mincov are blastn-only; every minimap2 entry point rejects
    them instead of silently ignoring them (reads thresholds have their own
    flags). ``mode`` names the invocation in the frozen message."""
    if minid != 80.0 or mincov != 80.0:
        usage_fail(
            f"--minid/--mincov apply to blastn only; use --min-identity/--min-breadth with {mode}"
        )


def _validate_reads_usage(
    output_format: OutputFormat | None,
    min_breadth: float,
    min_identity: float,
    min_mapq: int,
    threads: int,
) -> None:
    """Usage gates shared by both minimap2 entry points, in the frozen order
    (format first, then thresholds), before any file is touched."""
    if output_format is OutputFormat.tsv or output_format is OutputFormat.csv:
        usage_fail("--format tsv|csv is not available in reads mode (use json or md)")
    if not 0.0 <= min_breadth <= 100.0:
        usage_fail(f"--min-breadth must be in [0, 100]: got {min_breadth}")
    if not 0.0 <= min_identity <= 100.0:
        usage_fail(f"--min-identity must be in [0, 100]: got {min_identity}")
    if min_mapq < 0:
        usage_fail(f"--min-mapq must be >= 0: got {min_mapq}")
    if threads < 1:
        usage_fail(f"--threads must be >= 1: got {threads}")


def _screen_lanes(
    lanes: list[tuple[Path, Path | None]],
    db_name: str,
    datadir: Path | None,
    read_type: ReadTypeEnum | None,
    min_breadth: float,
    min_identity: float,
    min_mapq: int,
    threads: int,
    output_format: OutputFormat | None,
    quiet: bool,
    debug: bool,
) -> str:
    """Minimap2 engine core shared by both entry points: preset resolution,
    screening, rendering; json is the default format (SPEC.md §10). Either
    reads/2 threshold on selects the gapit.reads/2 document; both off keep
    gapit.reads/1 byte-identical. Returns the rendered output for the caller
    to echo."""
    resolved = _resolve_read_preset(lanes, read_type, quiet)
    database = find_database(config.resolve_datadir(datadir), db_name)
    read_files = [r1_path for r1_path, _ in lanes] + [
        r2_path for _, r2_path in lanes if r2_path is not None
    ]
    read_list = ", ".join(str(path) for path in read_files)
    if not quiet:
        typer.echo(f"Screening reads: {read_list}", err=True)
    report = screen_reads(
        lanes,
        database,
        read_type=resolved.value,
        min_breadth=min_breadth,
        threads=threads,
        debug=debug,
        min_identity=min_identity,
        min_mapq=min_mapq,
    )
    present = sum(1 for gene in report.genes if gene.present)
    if not quiet:
        typer.echo(f"Detected {present} present genes in {read_list}", err=True)
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
    if reads2:
        output = (
            render_reads2_markdown([report], params, now=now)
            if output_format is OutputFormat.md
            else render_reads2_json([report], params, now=now)
        )
    else:
        output = (
            render_reads_markdown([report], params, now=now)
            if output_format is OutputFormat.md
            else render_reads_json([report], params, now=now)
        )
    return output


def run_screen_reads(
    r1: list[Path],
    r2: list[Path] | None,
    db_name: str,
    datadir: Path | None,
    read_type: ReadTypeEnum | None,
    min_breadth: float,
    min_identity: float,
    min_mapq: int,
    threads: int,
    output_format: OutputFormat | None,
    quiet: bool,
    debug: bool = False,
    aligner: AlignerEnum | None = None,
    minid: float = 80.0,
    mincov: float = 80.0,
) -> str:
    """Screen FASTQ reads or assembly FASTA given as already-split per-lane
    --r1/--r2 file lists (the CLI owns the comma-splitting; MCP passes arrays
    natively, so commas in filenames survive). Per-lane minimap2, sample-level
    union; json is the default format (SPEC.md §10). A nonzero
    --min-identity/--min-mapq turns on gapit.reads/2 alignment filtering.
    Returns the rendered output."""
    if aligner is AlignerEnum.blastn:
        usage_fail("--aligner blastn is not available for --r1/--r2 reads input")
    _validate_reads_usage(output_format, min_breadth, min_identity, min_mapq, threads)
    _reject_blastn_thresholds(minid, mincov, "--r1/--r2")
    lanes = _pair_read_lanes(r1, r2)
    for path in [r1_path for r1_path, _ in lanes] + [
        r2_path for _, r2_path in lanes if r2_path is not None
    ]:
        ensure_input_file(path, "reads file")
    return _screen_lanes(
        lanes,
        db_name,
        datadir,
        read_type,
        min_breadth,
        min_identity,
        min_mapq,
        threads,
        output_format,
        quiet,
        debug,
    )


def run_screen_assemblies(
    files: list[Path] | None,
    fofn: Path | None,
    db_name: str,
    datadir: Path | None,
    read_type: ReadTypeEnum | None,
    min_breadth: float,
    min_identity: float,
    min_mapq: int,
    threads: int,
    jobs: int,
    noheader: bool,
    nopath: bool,
    output_format: OutputFormat | None,
    quiet: bool,
    debug: bool = False,
    minid: float = 80.0,
    mincov: float = 80.0,
) -> str:
    """Screen positional assembly FASTA file(s) with the minimap2 engine
    (--aligner minimap2): every input must be FASTA(.gz) content — FASTQ
    content is a usage error, undetectable content keeps the typed input
    error. Preset resolution and output follow the reads contract (SPEC §10);
    the blastn-engine-only flags --fofn/--jobs/--noheader/--nopath and the
    blastn thresholds --minid/--mincov are rejected here instead of silently
    ignored. A nonzero --min-identity/--min-mapq turns on gapit.reads/2
    filtering. Returns the rendered output."""
    _validate_reads_usage(output_format, min_breadth, min_identity, min_mapq, threads)
    if fofn is not None:
        usage_fail("--fofn is not available with --aligner minimap2")
    if jobs != 1:
        usage_fail("--jobs is not available with --aligner minimap2")
    if noheader:
        usage_fail("--noheader is not available with --aligner minimap2")
    if nopath:
        usage_fail("--nopath is not available with --aligner minimap2")
    _reject_blastn_thresholds(minid, mincov, "--aligner minimap2")
    if not files:
        usage_fail("no input files given (positional FILEs)")
    for path in files:
        ensure_input_file(path)
        if detect_read_kind(path) is not ReadFileKind.fasta:
            usage_fail("minimap2 engine requires FASTA assemblies")
    return _screen_lanes(
        [(path, None) for path in files],
        db_name,
        datadir,
        read_type,
        min_breadth,
        min_identity,
        min_mapq,
        threads,
        output_format,
        quiet,
        debug,
    )
