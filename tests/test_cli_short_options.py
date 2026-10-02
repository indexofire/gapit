"""CLI tests for short-form option aliases (--output/-o style pairs).

Two layers over the same contract:
- structural: the built command tree must give every option a short form;
  shorts are single letters, unique within their command, never -h;
- behavioral: commands invoked with shorts behave exactly like the long
  spelling — TSV/CSV/receipts byte-compared, JSON compared after dropping
  the wall-clock tool.created_at metadata field.
"""

import hashlib
import json
import re
import shutil
from collections.abc import Iterator
from pathlib import Path

import pytest
from typer._click import Command  # typer's vendored click: what CliRunner parses with
from typer.core import TyperGroup
from typer.main import get_command
from typer.testing import CliRunner, Result

from gapit.cli import app
from gapit.db import make_blast_db
from gapit.records import Manifest, Record, write_manifest, write_records

DATA = Path(__file__).parent / "data"
CONTIGS = DATA / "contigs"
READS = DATA / "reads"
CLUSTER = DATA / "cluster"
SUMMARY = DATA / "summary"
SEQ = "ACGTAG" * 10

runner = CliRunner()
ANSI_STYLE = re.compile(r"\x1b\[[0-9;]*m")


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    """Fresh datadir with a built tinyamr BLAST index (test_cli_screen.py)."""
    target = tmp_path / "datadir"
    shutil.copytree(DATA / "db", target)
    make_blast_db(target / "tinyamr" / "sequences", "tinyamr")
    return target


@pytest.fixture()
def records_datadir(tmp_path: Path) -> Path:
    """One hand-written installed database (test_cli_db_query.py convention)."""
    db_dir = tmp_path / "records"
    db_dir.mkdir()
    write_records(
        (Record(db="myamr", gene="tetX", sequence=SEQ, product="demo"),),
        db_dir / "records.jsonl",
    )
    write_manifest(
        Manifest(
            name="myamr",
            source_urls=(),
            fetched_at="2020-01-01T00:00:00+00:00",
            sha256="0" * 64,
            n_records=1,
            dbtype="nucl",
        ),
        db_dir / "gapit-manifest.json",
    )
    return tmp_path


def _collapse(result: Result) -> str:
    """ANSI-stripped, box-border-free, whitespace-collapsed help text."""
    text = ANSI_STYLE.sub("", result.stdout)
    for box_char in "│├└┼╭╮╰╯─":
        text = text.replace(box_char, " ")
    return " ".join(text.split())


def _json_without_clock(raw: str) -> str:
    """Stable-comparable JSON: the metadata wall-clock field is dropped
    (top-level `created_at`; second-boundary straddling otherwise flakes)."""
    document = json.loads(raw)
    document.pop("created_at", None)
    document.get("tool", {}).pop("created_at", None)
    return json.dumps(document, sort_keys=True)


def _assert_same_run(short: Result, long: Result) -> None:
    assert short.exit_code == long.exit_code, (short.stderr, long.stderr)
    assert short.stdout == long.stdout


def _walk(cmd: Command) -> Iterator[Command]:
    yield cmd
    if isinstance(cmd, TyperGroup):
        for sub in cmd.commands.values():
            yield from _walk(sub)


def test_every_option_has_a_unique_single_letter_short() -> None:
    """Given the built command tree, When every command's options are
    inspected, Then each long form carries a one-letter short, unique within
    its command and never -h (click's --help alias). The only exceptions are
    the typer-managed completion flags, which stay long-only by directive."""
    typer_managed = {"--install-completion", "--show-completion"}
    seen = 0
    for cmd in _walk(get_command(app)):
        taken: dict[str, str] = {}
        for param in cmd.params:
            longs = [o for o in param.opts if o.startswith("--")]
            if not longs or longs[0] in typer_managed:
                continue  # positional argument, or a typer-managed completion flag
            shorts = [o for o in param.opts if o.startswith("-") and not o.startswith("--")]
            assert shorts, f"{cmd.name}: {longs[0]} lacks a short form"
            (short,) = shorts
            assert len(short) == 2, f"{cmd.name}: {longs[0]} short is not a single letter"
            assert short != "-h", f"{cmd.name}: {longs[0]} would collide with --help/-h"
            assert short not in taken, (
                f"{cmd.name}: {short} is claimed by both {taken[short]} and {longs[0]}"
            )
            taken[short] = longs[0]
            seen += 1
    assert seen >= 40  # the walk really visited the whole command tree


def test_screen_help_renders_paired_shorts() -> None:
    """Given `gapit screen -h`, When rendered, Then long and short forms
    appear side by side and the negative --no-merge-fragments half stays
    long-only."""
    result = runner.invoke(app, ["screen", "-h"], env={"COLUMNS": "200"})
    assert result.exit_code == 0
    collapsed = _collapse(result)
    for pair in (
        "--db -d",
        "--datadir -D",
        "--format -f",
        "--output -o",
        "--threads -t",
        "--quiet -q",
        "--r1 -1",
        "--fofn -F",
        "--merge-fragments -m",
    ):
        assert pair in collapsed, pair
    assert "--no-merge-fragments" in collapsed


def test_screen_short_flags_match_long_form(datadir: Path) -> None:
    """Given a tsv screen run spelled with shorts, When compared to the long
    spelling, Then stdout is byte-identical."""
    contig = str(CONTIGS / "sort.fa")
    _assert_same_run(
        runner.invoke(app, ["screen", contig, "-d", "tinyamr", "-D", str(datadir), "-p"]),
        runner.invoke(
            app, ["screen", contig, "--db", "tinyamr", "--datadir", str(datadir), "--nopath"]
        ),
    )


def test_screen_short_format_json_matches_long_form(datadir: Path) -> None:
    """Given -f json vs --format json, When compared after dropping the
    metadata clock, Then the documents are identical."""
    contig = str(CONTIGS / "sort.fa")
    short = runner.invoke(
        app, ["screen", contig, "-d", "tinyamr", "-D", str(datadir), "-f", "json"]
    )
    long = runner.invoke(
        app, ["screen", contig, "--db", "tinyamr", "--datadir", str(datadir), "--format", "json"]
    )
    assert short.exit_code == long.exit_code == 0, (short.stderr, long.stderr)
    assert _json_without_clock(short.stdout) == _json_without_clock(long.stdout)


def test_screen_short_output_writes_file(datadir: Path, tmp_path: Path) -> None:
    """Given -o PATH, When run, Then the file holds the same bytes as a
    long-form --output run and as a plain stdout run, and stdout stays empty."""
    contig = str(CONTIGS / "sort.fa")
    short_path, long_path = tmp_path / "short.tsv", tmp_path / "long.tsv"
    short = runner.invoke(
        app, ["screen", contig, "-d", "tinyamr", "-D", str(datadir), "-o", str(short_path)]
    )
    long = runner.invoke(
        app,
        [
            "screen",
            contig,
            "--db",
            "tinyamr",
            "--datadir",
            str(datadir),
            "--output",
            str(long_path),
        ],
    )
    plain = runner.invoke(app, ["screen", contig, "--db", "tinyamr", "--datadir", str(datadir)])
    assert short.exit_code == long.exit_code == 0, short.stderr
    assert short.stdout == ""
    assert short_path.read_bytes() == long_path.read_bytes() == plain.stdout.encode()


def test_screen_short_reads_r1_r2_match_long_form(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given reads mode spelled with -1/-2, When compared to --r1/--r2, Then
    the gapit.reads documents match (clock field dropped)."""
    datadir = tmp_path / "datadir"
    shutil.copytree(DATA / "reads_db", datadir)
    monkeypatch.chdir(READS)
    short = runner.invoke(
        app,
        [
            "screen",
            "-1",
            "tetx_pe1_R1.fq",
            "-2",
            "tetx_pe1_R2.fq",
            "-d",
            "tinyreads",
            "-D",
            str(datadir),
        ],
    )
    long = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "tetx_pe1_R1.fq",
            "--r2",
            "tetx_pe1_R2.fq",
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
        ],
    )
    assert short.exit_code == long.exit_code == 0, (short.stderr, long.stderr)
    assert _json_without_clock(short.stdout) == _json_without_clock(long.stdout)


def test_screen_short_quiet_matches_long_form(datadir: Path) -> None:
    """Given -q vs --quiet, When compared, Then both silence stderr
    identically while stdout stays byte-identical."""
    contig = str(CONTIGS / "sort.fa")
    short = runner.invoke(app, ["screen", contig, "-d", "tinyamr", "-D", str(datadir), "-q"])
    long = runner.invoke(
        app, ["screen", contig, "--db", "tinyamr", "--datadir", str(datadir), "--quiet"]
    )
    _assert_same_run(short, long)
    assert short.stderr == ""
    plain = runner.invoke(app, ["screen", contig, "--db", "tinyamr", "--datadir", str(datadir)])
    assert "Processing:" in plain.stderr


def test_db_build_short_kind_typing_match_long_form(tmp_path: Path) -> None:
    """Given a cluster build spelled with -k/-T (and -D), When compared to
    the long spelling, Then the sequence, feature, and typing artifacts are
    byte-identical (the manifest alone differs: wall-clock fetched_at)."""
    for spelling, dd, dk, dt in (
        ("short", "-D", "-k", "-T"),
        ("long", "--datadir", "--kind", "--typing"),
    ):
        datadir = tmp_path / spelling
        datadir.mkdir()
        result = runner.invoke(
            app,
            [
                "db",
                "build",
                "myclu",
                str(CLUSTER / "bakta_style.gbk"),
                dd,
                str(datadir),
                dk,
                "cluster",
                dt,
                str(CLUSTER / "typing_valid.json"),
            ],
        )
        assert result.exit_code == 0, result.stderr
    for artifact in ("sequences", "features.json", "typing.json"):
        assert (tmp_path / "short" / "myclu" / artifact).read_bytes() == (
            tmp_path / "long" / "myclu" / artifact
        ).read_bytes(), artifact


def test_db_list_short_json_matches_long_form(records_datadir: Path) -> None:
    """Given db list -J vs --json, When compared, Then stdout is byte-identical."""
    _assert_same_run(
        runner.invoke(app, ["db", "list", "-D", str(records_datadir), "-J"]),
        runner.invoke(app, ["db", "list", "--datadir", str(records_datadir), "--json"]),
    )


def test_db_search_short_flags_match_long_form(records_datadir: Path) -> None:
    """Given db search spelled entirely with shorts, When compared to the
    long spelling, Then stdout is byte-identical."""
    _assert_same_run(
        runner.invoke(
            app,
            [
                "db",
                "search",
                "tetx",
                "-D",
                str(records_datadir),
                "-d",
                "myamr",
                "-f",
                "gene",
                "-e",
                "-l",
                "5",
            ],
        ),
        runner.invoke(
            app,
            [
                "db",
                "search",
                "tetx",
                "--datadir",
                str(records_datadir),
                "--db",
                "myamr",
                "--field",
                "gene",
                "--exact",
                "--limit",
                "5",
            ],
        ),
    )


def test_db_outdated_short_flags_match_long_form(records_datadir: Path) -> None:
    """Given db outdated -d/-J vs --days/--json, When compared, Then stdout
    is byte-identical."""
    _assert_same_run(
        runner.invoke(app, ["db", "outdated", "-D", str(records_datadir), "-d", "30", "-J"]),
        runner.invoke(
            app,
            ["db", "outdated", "--datadir", str(records_datadir), "--days", "30", "--json"],
        ),
    )


def test_db_install_short_flags_match_long_form(tmp_path: Path) -> None:
    """Given db install -s/-o vs --sha256/--output, When both target the
    same destination, Then the receipts are byte-identical and the installed
    bytes match the source."""
    source = tmp_path / "source.bin"
    destination = tmp_path / "installed.bin"
    source.write_bytes(b"short-form alias bytes\n")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    _assert_same_run(
        runner.invoke(app, ["db", "install", str(source), "-s", digest, "-o", str(destination)]),
        runner.invoke(
            app,
            [
                "db",
                "install",
                str(source),
                "--sha256",
                digest,
                "--output",
                str(destination),
            ],
        ),
    )
    assert destination.read_bytes() == source.read_bytes()


def test_root_short_version_and_json_match_long_form() -> None:
    """Given -V (and -V -J) vs --version (--version --json), When compared,
    Then stdout is byte-identical."""
    _assert_same_run(runner.invoke(app, ["-V"]), runner.invoke(app, ["--version"]))
    _assert_same_run(runner.invoke(app, ["-V", "-J"]), runner.invoke(app, ["--version", "--json"]))


def test_summary_short_format_matches_long_form(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given summary -f csv vs --format csv, When compared, Then stdout is
    byte-identical."""
    monkeypatch.chdir(SUMMARY)
    _assert_same_run(
        runner.invoke(app, ["summary", "-f", "csv", "sample_a.csv", "sample_b.csv"]),
        runner.invoke(app, ["summary", "--format", "csv", "sample_a.csv", "sample_b.csv"]),
    )
