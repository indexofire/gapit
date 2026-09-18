"""CLI tests for shell completions (typer add_completion on the root app)."""

from typer.testing import CliRunner

from gapit.cli import app

runner = CliRunner()


def test_root_help_lists_completion_options() -> None:
    """Given --help on the root app, When inspected, Then both typer
    completion options are offered."""
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "--install-completion" in result.stdout
    assert "--show-completion" in result.stdout


def test_show_completion_emits_shell_script() -> None:
    """Given --show-completion under a pinned SHELL, When run, Then a bash
    completion script referencing the gapit completion env var is printed
    to stdout and nothing else leaks into other runs."""
    result = runner.invoke(app, ["--show-completion"], env={"SHELL": "/bin/bash"})
    assert result.exit_code == 0
    assert "_GAPIT_COMPLETE" in result.stdout
    assert "complete -o default -F _gapit_completion gapit" in result.stdout


def test_bash_completion_lists_subcommands() -> None:
    """Given the bash completion protocol invocation (COMP_WORDS at the
    subcommand position), When run, Then the known subcommand names are
    emitted on stdout."""
    result = runner.invoke(
        app,
        [],
        env={
            "_GAPIT_COMPLETE": "complete_bash",
            "COMP_WORDS": "gapit ",
            "COMP_CWORD": "1",
        },
    )
    assert result.exit_code == 0
    completions = result.stdout.splitlines()
    assert "screen" in completions
    assert "setupdb" in completions
    assert "summary" in completions


def test_completion_output_only_when_requested() -> None:
    """Given a plain --version run, When run, Then stdout carries only the
    version line (completion machinery stays silent unless invoked)."""
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert "_gapit_completion" not in result.stdout
