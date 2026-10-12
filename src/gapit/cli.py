"""gapit command-line interface (click entrypoint)."""

import json
from pathlib import Path

import click

from gapit import __version__, config, db
from gapit.bundled import bundled_names, materialize_bundled
from gapit.cmd_db import register_db_command
from gapit.cmd_screen import register_screen_command
from gapit.cmd_summary import register_summary_command
from gapit.cmd_typing import register_typing_command
from gapit.dispatch import datadir_option, dispatch
from gapit.errors import usage_fail
from gapit.formats.json import VersionDocument
from gapit.formats.schemas import SCHEMA_MODELS
from gapit.mcp import register_mcp_command


@click.group(
    name="gapit",
    help="Mass screening of contigs and reads for known genes across reference databases.",
    no_args_is_help=True,
    invoke_without_command=True,
    # -h alias for --help; click Context propagates this to all subcommands/sub-apps.
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.option(
    "--version",
    "-V",
    "show_version",
    is_flag=True,
    default=None,
    help="Show version and exit.",
)
@click.option(
    "--json",
    "-J",
    "as_json",
    is_flag=True,
    help="With --version: emit gapit.version/1 JSON.",
)
@click.pass_context
def main(ctx: click.Context, show_version: bool | None, as_json: bool) -> None:
    """Mass screening of contigs and reads for known genes across reference databases."""
    if show_version:
        if as_json:
            click.echo(VersionDocument(version=__version__).model_dump_json(by_alias=True))
        else:
            click.echo(f"gapit {__version__}")
        ctx.exit(0)
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())
        ctx.exit(0)


app = main


def _setupdb(datadir: Path | None, debug: bool) -> None:
    root = config.resolve_datadir(datadir)
    for bundled_name in bundled_names():
        if not (root / bundled_name / "gapit-manifest.json").is_file():
            materialize_bundled(bundled_name, root, quiet=False)
    infos = db.list_databases(root, setupdb=True, debug=debug)
    for info in infos:
        click.echo(
            f"Indexed {info.name} ({info.n_sequences} sequences, {info.dbtype})",
            err=True,
        )


@app.command("setupdb")
@datadir_option()
@click.option("--debug", "-v", is_flag=True, help="Echo external command lines to stderr.")
def setupdb(datadir: Path | None, debug: bool) -> None:
    """Materialize bundled databases, then build BLAST indices for everything under the datadir."""
    dispatch(lambda: _setupdb(datadir, debug))


register_screen_command(app)
register_typing_command(app)
register_summary_command(app)
register_db_command(app)
register_mcp_command(app)


@app.command("schema", no_args_is_help=True)
@click.argument(
    "name",
    help=(
        "Document to introspect: report, reads, reads2, cluster, cluster2,"
        " summary, error, version, features, typing, typing_result, or floors."
    ),
)
def schema(name: str) -> None:
    """Print the JSON Schema of a gapit output document."""

    def run() -> None:
        model = SCHEMA_MODELS.get(name)
        if model is None:
            usage_fail(f"unknown schema name: {name} (choose from: {', '.join(SCHEMA_MODELS)})")
        click.echo(json.dumps(model.model_json_schema(by_alias=True), indent=2))

    dispatch(run)
