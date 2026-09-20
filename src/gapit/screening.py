"""The blastn contig-screening use-case plus helpers shared by both engines.

The minimap2 engine use-cases (--r1/--r2 reads, --aligner minimap2 assemblies)
live in screening_reads.py; this module keeps the abricate-parity contig
pipeline and the shared OutputFormat / AlignerEnum / database lookup.
"""

import enum
import os
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from typing import NoReturn

import typer

from gapit import config, db
from gapit.blast import ensure_blast, screen_file
from gapit.errors import DatabaseError, InputError, UsageError
from gapit.formats.json import render_json
from gapit.formats.md import render_markdown
from gapit.formats.tsv import format_tsv
from gapit.report import Report, ScreeningParams


class OutputFormat(enum.Enum):
    """Screen output formats."""

    tsv = "tsv"
    csv = "csv"
    json = "json"
    md = "md"


class AlignerEnum(enum.Enum):
    """Alignment engines for screen (SPEC.md §1/§10)."""

    blastn = "blastn"
    minimap2 = "minimap2"


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


def find_database(datadir: Path, name: str) -> db.Database:
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
    noheader: bool,
    nopath: bool,
    debug: bool,
    output_format: OutputFormat,
) -> None:
    """Screen each input file in order; buffer reports; print once at the end.

    The per-run gates (blastn presence via ``ensure_blast``, dbtype via one
    ``blastdbcmd -info``) fire once up front, so MISSING_DEPENDENCY and
    DATABASE_NOT_INDEXED errors precede any "Processing:" stderr lines.
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
    params = ScreeningParams(db=db_name, minid=minid, mincov=mincov, threads=threads)
    database = find_database(config.resolve_datadir(datadir), db_name)
    ensure_blast()
    dbtype = db.blast_db_info(database.sequences_path).dbtype
    cpu_count = os.cpu_count()
    if not quiet and cpu_count is not None and jobs * threads > cpu_count:
        typer.echo(f"--jobs {jobs} --threads {threads} oversubscribes {cpu_count} cpus", err=True)
    reports: list[Report] = []
    if jobs == 1:
        for path in inputs:
            if not quiet:
                typer.echo(f"Processing: {path}", err=True)
            report = screen_file(path, database, params, dbtype=dbtype, debug=debug)
            if not quiet:
                typer.echo(f"Found {len(report.hits)} genes in {path}", err=True)
            reports.append(report)
    else:

        def screen_one(path: Path) -> Report:
            if not quiet:
                typer.echo(f"Processing: {path}", err=True)
            report = screen_file(path, database, params, dbtype=dbtype, debug=debug)
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
        as_csv = output_format is OutputFormat.csv
        typer.echo(format_tsv(reports, csv=as_csv, noheader=noheader, nopath=nopath), nl=False)
