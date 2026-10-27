"""The `gapit db info` command: one database's provenance and guidance.

Read-only: it reports what is installed, what ships in the wheel, or what the
provider catalog declares — including the literature the content came from and
the caveats its interpretation depends on. ``--json`` emits the
``gapit.dbinfo/1`` document; the default is a stable ``FIELD<TAB>VALUE``
table so both humans and agents can read it without a TTY fork.
"""

from pathlib import Path

import click

from gapit.db_info_ops import collect_db_info, render_db_info
from gapit.dispatch import datadir_option, dispatch


@click.argument("name", help="Database name (see: gapit db list).")
@datadir_option()
@click.option(
    "--json", "-J", "as_json", is_flag=True, help="Print the gapit.dbinfo/1 JSON document."
)
def db_info_command(name: str, datadir: Path | None, as_json: bool) -> None:
    """Show one database's provenance, references, usage, and caveats."""

    def run() -> None:
        document = collect_db_info(name, datadir)
        if as_json:
            click.echo(document.model_dump_json(by_alias=True, exclude_none=True))
            return
        click.echo(render_db_info(document), nl=False)

    dispatch(run)
