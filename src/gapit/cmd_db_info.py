"""The `gapit db info` command: one database's provenance and guidance.

Read-only: it reports what is installed, what ships in the wheel, or what the
provider catalog declares — including the literature the content came from and
the caveats its interpretation depends on.

Three renderings, picked like the rest of the CLI (``db list`` precedent):

- **manual** (default on a terminal) — a plain man-style page: NAME, SYNOPSIS,
  DESCRIPTION, STATUS, PROVENANCE, REFERENCES, CAVEATS, USAGE, SEE ALSO;
- **tsv** (default when stdout is piped, or forced with ``--tsv``) — the stable
  two-column ``FIELD<TAB>VALUE`` stream agents already parse;
- **json** (``--json``) — the ``gapit.dbinfo/1`` document.
"""

import sys
from pathlib import Path

import click

from gapit.db_info_manual import render_db_manual
from gapit.db_info_ops import collect_db_info, render_db_info
from gapit.dispatch import datadir_option, dispatch


def _stdout_is_tty() -> bool:
    """Whether the terminal (not a pipe) is reading — the manual/table fork."""
    return sys.stdout.isatty()


@click.option(
    "--db",
    "-d",
    "db_name",
    required=True,
    help="Database name to describe (see: gapit db list).",
)
@datadir_option()
@click.option(
    "--json", "-J", "as_json", is_flag=True, help="Print the gapit.dbinfo/1 JSON document."
)
@click.option(
    "--tsv", "-t", "as_tsv", is_flag=True, help="Force the FIELD<TAB>VALUE table (pipe default)."
)
def db_info_command(db_name: str, datadir: Path | None, as_json: bool, as_tsv: bool) -> None:
    """Show one database's provenance, references, usage, and caveats.

    On a terminal this prints a man-style page; piped output stays the stable
    FIELD<TAB>VALUE table (forceable with --tsv), and --json emits
    gapit.dbinfo/1.
    """

    def run() -> None:
        document = collect_db_info(db_name, datadir)
        if as_json:
            click.echo(document.model_dump_json(by_alias=True, exclude_none=True))
            return
        if as_tsv or not _stdout_is_tty():
            click.echo(render_db_info(document), nl=False)
            return
        click.echo(render_db_manual(document), nl=False)

    dispatch(run)
