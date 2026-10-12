"""The `gapit db build` command: click shell over the use-case
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

import click

from gapit.db_build_ops import Dbtype, Kind, perform_build
from gapit.dispatch import datadir_option, dispatch, quiet_option


@click.argument("name", help="Target database name (created under the datadir).")
@click.argument(
    "fasta",
    type=click.Path(path_type=Path),
    help=(
        "Input FASTA (gene db) or GBK/GFF3 file (cluster db), detected by suffix"
        " (.gz/.bz2 accepted)."
    ),
)
@click.option(
    "--tsv",
    "-t",
    type=click.Path(path_type=Path),
    help="Metadata TSV: header row with gene/accession/function columns.",
)
@click.option(
    "--dbtype",
    "-d",
    type=click.Choice(("nucl", "prot")),
    help="Force nucl or prot (default: auto-detect).",
)
@click.option(
    "--kind",
    "-k",
    type=click.Choice(("gene", "cluster")),
    help="Force gene or cluster (default: auto-detect by suffix; must agree with it).",
)
@click.option(
    "--typing",
    "-T",
    type=click.Path(path_type=Path),
    help=(
        "gapit.typing/1 or /2 scoring spec: validated, copied into the db,"
        " recorded in its manifest."
    ),
)
@click.option(
    "--floors",
    "-F",  # uppercase: -f is taken by --force.
    type=click.Path(path_type=Path),
    help=(
        "gapit.floors/1 per-gene identity floors for reads-mode presence:"
        " validated, copied into the db as floors.json (gene builds only)."
    ),
)
@datadir_option()
@click.option(
    "--description",
    "-e",  # second letter: -d/-D are taken by --dbtype/--datadir.
    default="",
    help="Default product for records whose FASTA header has no description text.",
)
@click.option("--force", "-f", is_flag=True, help="Overwrite the database if it already exists.")
@quiet_option()
def db_build_command(
    name: str,
    fasta: Path,
    tsv: Path | None,
    dbtype: Dbtype | None,
    kind: Kind | None,
    typing: Path | None,
    floors: Path | None,
    datadir: Path | None,
    description: str,
    force: bool,
    quiet: bool,
) -> None:
    """Build a custom gapit-native database from a FASTA/GBK/GFF3 input."""

    def run() -> None:
        def warn(message: str) -> None:
            if not quiet:
                click.echo(f"WARNING: {message}", err=True)

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
        click.echo(receipt.model_dump_json(exclude_none=True))

    dispatch(run)
