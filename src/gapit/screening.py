"""The screen use-case: input resolution, validation, and report rendering."""

import enum
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from typing import NoReturn

import typer

from gapit import config, db
from gapit.blast import screen_file
from gapit.errors import DatabaseError, InputError, UsageError
from gapit.formats.json import render_json, render_reads_json
from gapit.formats.md import render_markdown, render_reads_markdown
from gapit.formats.tsv import format_tsv
from gapit.reads import ReadsParams, ReadTypeEnum, screen_reads
from gapit.report import Report, ScreeningParams


class OutputFormat(enum.Enum):
    """Screen output formats."""

    tsv = "tsv"
    csv = "csv"
    json = "json"
    md = "md"


def usage_fail(message: str) -> NoReturn:
    """Raise a usage error (gapit.error/1 envelope, exit 2)."""
    raise UsageError(message, code="USAGE_ERROR")


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
        if not path.is_file():
            raise InputError(
                f"input file not found or unreadable: {path}",
                code="INPUT_NOT_FOUND",
                context={"file": str(path)},
            )
    return inputs


def _find_database(datadir: Path, name: str) -> db.Database:
    """Look up a database by name under the datadir; unknown names list what exists."""
    databases = db.discover_databases(datadir)
    for database in databases:
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
    csv_flag: bool,
    noheader: bool,
    nopath: bool,
    debug: bool,
    output_format: OutputFormat,
) -> None:
    """Screen each input file in order; buffer reports; print once at the end."""
    if not 0.0 < minid <= 100.0:
        usage_fail(f"--minid must be in (0, 100]: got {minid}")
    if not 0.0 <= mincov <= 100.0:
        usage_fail(f"--mincov must be in [0, 100]: got {mincov}")
    if threads < 1:
        usage_fail(f"--threads must be >= 1: got {threads}")
    if jobs < 1:
        usage_fail(f"--jobs must be >= 1: got {jobs}")
    inputs = _resolve_inputs(files, fofn)
    params = ScreeningParams(db=db_name, minid=minid, mincov=mincov, threads=threads)
    database = _find_database(config.resolve_datadir(datadir), db_name)
    reports: list[Report] = []
    if jobs == 1:
        for path in inputs:
            if not quiet:
                typer.echo(f"Processing: {path}", err=True)
            report = screen_file(path, database, params, debug=debug)
            if not quiet:
                typer.echo(f"Found {len(report.hits)} genes in {path}", err=True)
            reports.append(report)
    else:

        def screen_one(path: Path) -> Report:
            if not quiet:
                typer.echo(f"Processing: {path}", err=True)
            report = screen_file(path, database, params, debug=debug)
            if not quiet:
                typer.echo(f"Found {len(report.hits)} genes in {path}", err=True)
            return report

        # Subprocess-bound work (GIL irrelevant); executor.map collects
        # positionally, so reports stay in input order (SPEC.md §4) and a
        # failing file raises at its position, like the sequential loop.
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            reports = list(executor.map(screen_one, inputs))
    if output_format is OutputFormat.json:
        typer.echo(render_json(reports, params, now=datetime.now(UTC)), nl=False)
    elif output_format is OutputFormat.md:
        typer.echo(render_markdown(reports, params, now=datetime.now(UTC)), nl=False)
    else:
        as_csv = csv_flag or output_format is OutputFormat.csv
        typer.echo(format_tsv(reports, csv=as_csv, noheader=noheader, nopath=nopath), nl=False)


def _parse_read_lanes(r1: str, r2: str | None) -> list[tuple[Path, Path | None]]:
    """Split comma-separated --r1/--r2 file lists into per-sample lanes."""
    r1_list = _split_read_list(r1, "--r1")
    r2_list = _split_read_list(r2, "--r2") if r2 is not None else None
    if r2_list is not None and len(r2_list) != len(r1_list):
        usage_fail(
            f"--r2 has {len(r2_list)} files but --r1 has {len(r1_list)} (lanes must pair up)"
        )
    if r2_list is None:
        return [(path, None) for path in r1_list]
    return list(zip(r1_list, r2_list, strict=True))


def _split_read_list(raw: str, flag: str) -> list[Path]:
    parts = [part.strip() for part in raw.split(",")]
    if any(not part for part in parts):
        usage_fail(f"{flag} contains an empty element: {raw!r}")
    return [Path(part) for part in parts]


def run_screen_reads(
    r1: str,
    r2: str | None,
    db_name: str,
    datadir: Path | None,
    read_type: ReadTypeEnum,
    min_breadth: float,
    threads: int,
    output_format: OutputFormat | None,
    quiet: bool,
    debug: bool = False,
) -> None:
    """Screen FASTQ reads (per-lane minimap2, sample-level union); json is the
    default format (SPEC.md §10)."""
    if output_format is OutputFormat.tsv or output_format is OutputFormat.csv:
        usage_fail("--format tsv|csv is not available in reads mode (use json or md)")
    if not 0.0 <= min_breadth <= 100.0:
        usage_fail(f"--min-breadth must be in [0, 100]: got {min_breadth}")
    if threads < 1:
        usage_fail(f"--threads must be >= 1: got {threads}")
    lanes = _parse_read_lanes(r1, r2)
    read_files = [r1_path for r1_path, _ in lanes] + [
        r2_path for _, r2_path in lanes if r2_path is not None
    ]
    for path in read_files:
        if not path.is_file():
            raise InputError(
                f"reads file not found or unreadable: {path}",
                code="INPUT_NOT_FOUND",
                context={"file": str(path)},
            )
    database = _find_database(config.resolve_datadir(datadir), db_name)
    read_list = ", ".join(str(path) for path in read_files)
    if not quiet:
        typer.echo(f"Screening reads: {read_list}", err=True)
    report = screen_reads(
        lanes,
        database,
        read_type=read_type.value,
        min_breadth=min_breadth,
        threads=threads,
        debug=debug,
    )
    present = sum(1 for gene in report.genes if gene.present)
    if not quiet:
        typer.echo(f"Detected {present} present genes in {read_list}", err=True)
    params = ReadsParams(
        db=db_name, read_type=read_type.value, min_breadth=min_breadth, threads=threads
    )
    now = datetime.now(UTC)
    if output_format is OutputFormat.md:
        typer.echo(render_reads_markdown([report], params, now=now), nl=False)
    else:
        typer.echo(render_reads_json([report], params, now=now), nl=False)
