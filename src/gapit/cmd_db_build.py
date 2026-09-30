"""The `gapit db build` command: typer shell over the use-case
(:mod:`gapit.db_build_ops`), which turns a user-supplied input into a
fully built gapit-native database. FASTA inputs build gene databases
(records.jsonl -> sequences + BLAST index + manifest, written last);
GBK/GFF inputs build cluster databases (locus FASTA + gapit.features/1
feature table + manifest ``kind: cluster``, plus an optional validated
gapit.typing/1 scoring spec via --typing).
"""

from pathlib import Path
from typing import Annotated

import typer

from gapit.db_build_ops import Dbtype, Kind, perform_build
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
                "Input FASTA (gene db) or GBK/GFF3 file (cluster db), detected by suffix"
                " (.gz/.bz2 accepted)."
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
    kind: Annotated[
        Kind | None,
        typer.Option(
            "--kind",
            help="Force gene or cluster (default: auto-detect by suffix; must agree with it).",
        ),
    ] = None,
    typing: Annotated[
        Path | None,
        typer.Option(
            "--typing", help="gapit.typing/1 scoring spec, validated and copied into a cluster db."
        ),
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
    """Build a custom gapit-native database from a FASTA/GBK/GFF3 input."""

    def run() -> None:
        def warn(message: str) -> None:
            if not quiet:
                typer.echo(f"WARNING: {message}", err=True)

        receipt = perform_build(
            name,
            fasta,
            tsv,
            dbtype,
            description,
            datadir,
            force,
            warn=warn,
            quiet=quiet,
            kind=kind,
            typing=typing,
        )
        typer.echo(receipt.model_dump_json())

    dispatch(run)
