"""The `gapit db rebuild` command: re-project a gene database from its records.

The curation workflow `records.jsonl` exists for: edit the truth (rename a gene,
fix a product, fold two alleles under one name), then rebuild — the sequences
FASTA, self-check, BLAST index, and manifest are regenerated in place, staged so
a broken edit leaves the installed database untouched.
"""

from pathlib import Path

import click

from gapit.db_rebuild_ops import perform_rebuild
from gapit.dispatch import datadir_option, dispatch, quiet_option


@click.argument("name", help="Installed gene database to rebuild (see: gapit db list).")
@datadir_option()
@quiet_option()
def db_rebuild_command(name: str, datadir: Path | None, quiet: bool) -> None:
    """Rebuild <datadir>/NAME's artifacts from its records.jsonl.

    Gene databases only: a cluster database's locus FASTA and features.json are
    generated from a source GBK/GFF, so rebuild it with `gapit db build` instead.
    """

    def run() -> None:
        receipt = perform_rebuild(name, datadir, quiet=quiet)
        click.echo(receipt.model_dump_json())

    dispatch(run)
