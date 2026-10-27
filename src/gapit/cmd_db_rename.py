"""The `gapit db rename` command: click shell over the use-case.

Renaming is a datadir maintenance operation, so it reuses the standard
acquisition guards: a plain target name, an installed source (its manifest is
the provenance the rebuild keeps), and no silent overwrite. One JSON receipt
goes to stdout; progress notes go to stderr.
"""

from pathlib import Path

import click

from gapit.db_rename_ops import perform_rename
from gapit.dispatch import datadir_option, dispatch, quiet_option


@click.argument("old", help="Installed database to rename (see: gapit db list).")
@click.argument("new", help="New database name (see the naming convention in docs/databases.md).")
@datadir_option()
@click.option("--force", "-f", is_flag=True, help="Replace the target if it already exists.")
@quiet_option()
def db_rename_command(old: str, new: str, datadir: Path | None, force: bool, quiet: bool) -> None:
    """Rename an installed database OLD to NEW.

    The database is rebuilt from its own records under the new name (gene
    databases re-header every record; cluster databases copy their locus
    artifacts and re-stamp the manifest), then the old directory is removed.
    A failure at any point leaves OLD untouched.
    """

    def run() -> None:
        receipt = perform_rename(old, new, datadir, force=force, quiet=quiet)
        click.echo(receipt.model_dump_json())

    dispatch(run)
