"""The `gapit db build` command: typer shell over the use-case
(:mod:`gapit.db_build_ops`), which turns a user-supplied input into a
fully built gapit-native database. FASTA inputs build gene databases
(records.jsonl -> sequences + BLAST index + manifest, written last);
GBK/GFF inputs build cluster databases (locus FASTA + gapit.features/1
feature table + manifest ``kind: cluster``). Both branches accept a
``--typing`` scoring spec (gapit.typing/1 or /2): validated against the
input's records, copied into the database as ``typing.json``, and recorded
in the manifest's ``typing_schema``.
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
            "--tsv", "-t", help="Metadata TSV: header row with gene/accession/function columns."
        ),
    ] = None,
    dbtype: Annotated[
        Dbtype | None,
        typer.Option("--dbtype", "-d", help="Force nucl or prot (default: auto-detect)."),
    ] = None,
    kind: Annotated[
        Kind | None,
        typer.Option(
            "--kind",
            "-k",
            help="Force gene or cluster (default: auto-detect by suffix; must agree with it).",
        ),
    ] = None,
    typing: Annotated[
        Path | None,
        typer.Option(
            "--typing",
            "-T",
            help=(
                "gapit.typing/1 or /2 scoring spec: validated, copied into the db,"
                " recorded in its manifest."
            ),
        ),
    ] = None,
    floors: Annotated[
        Path | None,
        typer.Option(
            # -F (uppercase): -f is taken by --force.
            "--floors",
            "-F",
            help=(
                "gapit.floors/1 per-gene identity floors for reads-mode presence:"
                " validated, copied into the db as floors.json (gene builds only)."
            ),
        ),
    ] = None,
    datadir: Annotated[
        Path | None,
        typer.Option(
            "--datadir",
            "-D",
            help="Database directory (default: $GAPIT_DATADIR, then ~/.local/share/gapit/db).",
        ),
    ] = None,
    description: Annotated[
        str,
        typer.Option(
            # -e (second letter): -d/-D are taken by --dbtype/--datadir.
            "--description",
            "-e",
            help="Default product for records whose FASTA header has no description text.",
        ),
    ] = "",
    force: Annotated[
        bool,
        typer.Option("--force", "-f", help="Overwrite the database if it already exists."),
    ] = False,
    quiet: Annotated[
        bool, typer.Option("--quiet", "-q", help="Silence stderr diagnostics.")
    ] = False,
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
            floors=floors,
        )
        typer.echo(receipt.model_dump_json(exclude_none=True))

    dispatch(run)
