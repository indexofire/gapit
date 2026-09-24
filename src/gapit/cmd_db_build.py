"""The `gapit db build` command: typer shell over the use-case
(:mod:`gapit.db_build_ops`), which turns a user-supplied FASTA into a
fully built gapit-native database (records.jsonl -> sequences + BLAST
index + manifest, written last).
"""

from pathlib import Path
from typing import Annotated

import typer

from gapit.db_build_ops import Dbtype, perform_build
from gapit.dispatch import dispatch


def db_build_command(
    name: Annotated[
        str,
        typer.Argument(help="Target database name (created under the datadir)."),
    ],
    fasta: Annotated[
        Path,
        typer.Argument(
            help=(
                "Input FASTA: plain, abricate ~~~, or gapit| headers, detected per"
                " record (.gz/.bz2 accepted)."
            ),
        ),
    ],
    tsv: Annotated[
        Path | None,
        typer.Option(
            "--tsv", help="Metadata TSV: header row with gene/accession/function columns."
        ),
    ] = None,
    dbtype: Annotated[
        Dbtype | None,
        typer.Option("--dbtype", help="Force nucl or prot (default: auto-detect)."),
    ] = None,
    datadir: Annotated[
        Path | None,
        typer.Option(
            "--datadir",
            help="Database directory (default: $GAPIT_DATADIR, then ~/.local/share/gapit/db).",
        ),
    ] = None,
    description: Annotated[
        str,
        typer.Option(
            "--description",
            help="Default product for records whose FASTA header has no description text.",
        ),
    ] = "",
    force: Annotated[
        bool,
        typer.Option("--force", help="Overwrite the database if it already exists."),
    ] = False,
    quiet: Annotated[bool, typer.Option("--quiet", help="Silence stderr diagnostics.")] = False,
) -> None:
    """Build a custom gapit-native database from a FASTA (+ optional TSV)."""

    def run() -> None:
        def warn(message: str) -> None:
            if not quiet:
                typer.echo(f"WARNING: {message}", err=True)

        receipt = perform_build(
            name, fasta, tsv, dbtype, description, datadir, force, warn=warn, quiet=quiet
        )
        typer.echo(receipt.model_dump_json())

    dispatch(run)
