"""CLI integration tests: uniform help-on-bare for input-requiring commands.

Commands that need input print their full help when invoked with no
arguments at all (typer ``no_args_is_help``) instead of erroring; commands
that are VALID with no arguments (`mcp`, `setupdb`, `db list`,
`db outdated`) must keep running bare. Two stdin-aware commands keep their
own tty-guarded bare-help instead: `typing` and `summary` read a piped
table when stdin is not a terminal (their suites pin the probe).
Flag-present-but-input-missing invocations keep their typed usage-error
envelopes (locked per-command in their own suites); bare `db fetch` never
reaches the use-case, so the former omitted-NAME default-set download
cannot start by accident.
"""

import re
from pathlib import Path

import pytest
from typer.testing import CliRunner, Result

from gapit.cli import app

runner = CliRunner()

# typer's vendored click signals no-args-is-help via NoArgsIsHelpError (a
# UsageError subclass), so help-on-bare exits 2 — the same code the bare
# root app has always used. Recorded reality, not an aspiration.
BARE_HELP_EXIT = 2

# Every full help render carries this row; click's missing-argument error
# (which also prints a "Usage:" line) does not. Whitespace is normalized
# before matching because rich wraps help rows at the terminal width.
HELP_SENTINEL = "Show this message and exit."


def _normalized(result: Result) -> str:
    # GITHUB_ACTIONS=true makes typer force terminal styling on help output,
    # splitting tokens with SGR runs — strip them before matching (the
    # test_cli_completion.py precedent) and pin COLUMNS at invocation.
    text = result.stdout + result.stderr
    text = re.compile(r"\x1b\[[0-9;]*m").sub("", text)
    for box_char in "│├└┼╭╮╰╯─":
        text = text.replace(box_char, " ")
    return " ".join(text.split())


BARE_COMMANDS = [
    pytest.param([], id="root"),
    pytest.param(["screen"], id="screen"),
    pytest.param(["schema"], id="schema"),
    pytest.param(["db"], id="db"),
    pytest.param(["db", "fetch"], id="db-fetch"),
    pytest.param(["db", "install"], id="db-install"),
    pytest.param(["db", "build"], id="db-build"),
    pytest.param(["db", "search"], id="db-search"),
]


@pytest.mark.parametrize(("argv",), BARE_COMMANDS)
def test_bare_invocation_prints_help(argv: list[str], monkeypatch: pytest.MonkeyPatch) -> None:
    """Given any input-requiring command, When invoked with no arguments at
    all, Then the full help prints with the no-args-is-help exit code and no
    gapit.error envelope (REGISTRY emptied so a bare `db fetch` that leaked
    past help could never start a real download)."""
    monkeypatch.setattr("gapit.db_ops.REGISTRY", {})
    result = runner.invoke(app, argv, env={"COLUMNS": "100"})
    combined = _normalized(result)
    assert result.exit_code == BARE_HELP_EXIT
    assert "Usage:" in combined
    assert HELP_SENTINEL in combined
    assert "gapit.error" not in combined


def test_screen_flag_without_files_keeps_usage_error() -> None:
    """Given `screen --db card` (a flag present, input missing), When run,
    Then the typed usage-error envelope still fires — no_args_is_help covers
    only the zero-argument case."""
    result = runner.invoke(app, ["screen", "--db", "card"])
    assert result.exit_code == 2
    assert "no input files given" in result.stderr
    assert HELP_SENTINEL not in _normalized(result)


def test_db_fetch_flag_without_name_is_usage_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given `db fetch --datadir D` with no NAME, When run, Then exit 2 with
    a USAGE_ERROR envelope asking for a NAME or `all` — flags present means
    help does not fire, and the default-set download is never implied."""
    monkeypatch.setattr("gapit.db_ops.REGISTRY", {})
    result = runner.invoke(app, ["db", "fetch", "--datadir", str(tmp_path)])
    assert result.exit_code == 2
    assert "gapit.error" in result.stderr
    assert "USAGE_ERROR" in result.stderr
    assert "all" in result.stderr
    assert HELP_SENTINEL not in _normalized(result)


def test_valid_no_arg_commands_do_not_print_help(tmp_path: Path) -> None:
    """Given commands that are valid with no arguments, When run bare, Then
    they execute (or emit their own documented envelopes) without ever
    printing help: db list, db outdated, setupdb."""
    empty = tmp_path / "datadir"
    empty.mkdir()

    listed = runner.invoke(app, ["db", "list", "--datadir", str(empty)])
    assert listed.exit_code == 0
    assert listed.stdout.startswith("NAME\tPROVIDER\tSTATUS\tDBTYPE\tDESCRIPTION")
    assert HELP_SENTINEL not in listed.stdout

    # Empty datadir: the documented DATADIR_EMPTY envelope (exit 4), not help.
    outdated = runner.invoke(app, ["db", "outdated", "--datadir", str(empty)])
    assert outdated.exit_code == 4
    assert "DATADIR_EMPTY" in outdated.stderr
    assert HELP_SENTINEL not in _normalized(outdated)

    setup = runner.invoke(app, ["setupdb", "--datadir", str(empty)])
    assert setup.exit_code == 0
    assert HELP_SENTINEL not in _normalized(setup)


def test_mcp_bare_serves_without_help() -> None:
    """Given the MCP stdio server, When `gapit mcp` runs bare with an
    immediately-closed stdin, Then it serves until EOF and exits 0 without
    printing help."""
    result = runner.invoke(app, ["mcp"], input="")
    assert result.exit_code == 0
    assert result.output == ""
