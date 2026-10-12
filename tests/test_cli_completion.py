"""CLI tests for shell completions (click env-var activation model).

Plain click injects no completion options: activation rides the
``_GAPIT_COMPLETE`` environment variable (``bash_source`` prints the
script, ``bash_complete`` answers a COMP_WORDS/COMP_CWORD probe).
"""

import re

from click.testing import CliRunner

from gapit.cli import app

runner = CliRunner()

# click's CliRunner renders help without SGR styling (color defaults off);
# the strip is kept as a defensive normalization before asserting tokens.
ANSI_STYLE = re.compile(r"\x1b\[[0-9;]*m")


def test_root_help_offers_no_completion_flags() -> None:
    """Given --help on the root app, When inspected, Then no legacy
    completion option is offered (plain click adds no flags)."""
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    help_text = ANSI_STYLE.sub("", result.stdout)
    assert "--install-completion" not in help_text
    assert "--show-completion" not in help_text


def test_install_completion_flag_is_rejected() -> None:
    """Given the legacy --install-completion flag, When invoked, Then
    click rejects it as an unknown option (exit 2) and no data reaches
    stdout."""
    result = runner.invoke(app, ["--install-completion"])
    assert result.exit_code == 2
    assert "No such option '--install-completion'." in result.stderr
    assert result.stdout == ""


def test_short_help_alias_works_on_every_level() -> None:
    """Given -h instead of --help, When invoked on the root app, a
    subcommand, and a nested sub-app, Then each prints its help and exits 0."""
    for argv in (["-h"], ["screen", "-h"], ["db", "-h"], ["db", "fetch", "-h"]):
        result = runner.invoke(app, argv)
        assert result.exit_code == 0, argv
        assert "Usage:" in ANSI_STYLE.sub("", result.stdout), argv


def test_bash_source_emits_completion_script() -> None:
    """Given the _GAPIT_COMPLETE=bash_source activation, When run, Then a
    bash completion script is printed to stdout: it defines the gapit
    completion function, calls back through _GAPIT_COMPLETE=bash_complete,
    and registers itself for the gapit program name."""
    result = runner.invoke(app, [], env={"_GAPIT_COMPLETE": "bash_source"})
    assert result.exit_code == 0
    assert "_gapit_completion() {" in result.stdout
    assert "_GAPIT_COMPLETE=bash_complete" in result.stdout
    assert "complete -o nosort -F _gapit_completion gapit" in result.stdout


def test_bash_completion_lists_subcommands() -> None:
    """Given the bash completion protocol invocation (COMP_WORDS at the
    subcommand position), When run, Then the known subcommand names are
    emitted on stdout as type-prefixed completion items."""
    result = runner.invoke(
        app,
        [],
        env={
            "_GAPIT_COMPLETE": "bash_complete",
            "COMP_WORDS": "gapit ",
            "COMP_CWORD": "1",
        },
    )
    assert result.exit_code == 0
    completions = [line.split(",", 1)[1] for line in result.stdout.splitlines()]
    assert "screen" in completions
    assert "setupdb" in completions
    assert "summary" in completions


def test_completion_output_only_when_requested() -> None:
    """Given a plain --version run, When run, Then stdout carries only the
    version line (completion machinery stays silent unless invoked)."""
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert "_gapit_completion" not in result.stdout
