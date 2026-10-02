"""The `gapit typing` command suite: input parsing, typed errors, formats,
and the two-stage equivalence lock.

Most parsing tests run over hand-written screen tables (the input contract
is the abricate 15-column shape); the equivalence class builds the markers
fixture db and proves the command's calls equal the inline engine's calls
over the screened Report — the TSV round-trip loses nothing. Goldens pin
the gapit.typing_result/1 JSON/MD renders at a pinned timestamp. The
fixture-level designation suites live in test_gene_typing*.py and
test_typing_schemes.py.
"""

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from gapit import cmd_typing
from gapit.cli import app
from gapit.cluster import GeneCall, load_typing
from gapit.db import discover_databases
from gapit.formats.typing_result import (
    TypingFileResult,
    render_typing_md,
    render_typing_result_json,
)
from gapit.typing_gene import evaluate_gene_calls
from gapit.typing_input import read_typing_input

DATA = Path(__file__).parent / "data" / "typing"
GOLDEN = Path(__file__).parent / "golden"
DB = "markers"
PINNED_NOW = datetime(2026, 10, 2, 12, 0, 0, tzinfo=UTC)

runner = CliRunner()

HEADER = (
    "#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t"
    "%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE"
)


def table(tmp_path: Path, name: str, rows: list[str], *, header: str = HEADER) -> Path:
    """A screen result table from data lines (tab-separated by default)."""
    path = tmp_path / name
    path.write_text("\n".join((header, *rows)) + "\n", encoding="utf-8")
    return path


def row(
    file: str = "s.fa",
    gene: str = "marker_a",
    database: str = DB,
    identity: str = "100.00",
    coverage: str = "100.00",
) -> str:
    return (
        f"{file}\tc1\t1\t90\t+\t{gene}\t1-90/90\t===============\t0/0\t"
        f"{coverage}\t{identity}\t{database}\t-\t-\t-"
    )


@pytest.fixture(scope="module")
def typed_db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    datadir = tmp_path_factory.mktemp("typing_cmd")
    result = runner.invoke(
        app,
        [
            "db",
            "build",
            DB,
            str(DATA / "markers.fa"),
            "--datadir",
            str(datadir),
            "--typing",
            str(DATA / "typing_markers_v2.json"),
        ],
    )
    assert result.exit_code == 0, result.stderr
    return datadir


@pytest.fixture(scope="module")
def untyped_db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    datadir = tmp_path_factory.mktemp("typing_cmd_untyped")
    result = runner.invoke(
        app, ["db", "build", DB, str(DATA / "markers.fa"), "--datadir", str(datadir)]
    )
    assert result.exit_code == 0, result.stderr
    return datadir


class TestParsing:
    def test_fold_keeps_the_best_row_per_gene(self, tmp_path: Path) -> None:
        """Given two rows for one gene with distinct identity/coverage, When
        the input parses, Then the GeneCall carries the lexicographically
        best (identity, coverage) row's values."""
        path = table(
            tmp_path,
            "fold.tsv",
            [
                row(identity="99.00", coverage="100.00"),
                row(identity="99.50", coverage="90.00"),
                row(identity="99.50", coverage="95.00"),
            ],
        )
        parsed = read_typing_input([path])
        assert parsed.database == DB
        (entry,) = parsed.files
        assert entry.file == "s.fa"
        assert entry.calls["marker_a"].identity_pct == 99.5
        assert entry.calls["marker_a"].coverage_pct == 95.0
        assert entry.calls["marker_a"].verdict == "present"

    def test_first_row_wins_an_exact_tie(self, tmp_path: Path) -> None:
        """Given two equal (identity, coverage) rows for one gene, When the
        input parses, Then the first row wins — the inline fold's strict
        comparison, so screen order stays deterministic."""
        path = table(
            tmp_path,
            "tie.tsv",
            [row(file="s.fa"), row(file="s.fa")],
        )
        (entry,) = read_typing_input([path]).files
        assert len(entry.calls) == 1

    def test_rows_merge_across_tables_by_file(self, tmp_path: Path) -> None:
        """Given two tables whose rows share a FILE value, When typed, Then
        the FILE folds once across both tables in first-appearance order."""
        first = table(tmp_path, "one.tsv", [row(gene="marker_a")])
        second = table(
            tmp_path, "two.tsv", [row(file="t.fa", gene="marker_b"), row(gene="marker_c")]
        )
        parsed = read_typing_input([first, second])
        assert [entry.file for entry in parsed.files] == ["s.fa", "t.fa"]
        assert set(parsed.files[0].calls) == {"marker_a", "marker_c"}

    def test_csv_tables_are_detected(self, tmp_path: Path) -> None:
        """Given a comma-separated screen table, When the input parses,
        Then the separator is auto-detected (summary.py rule) and the calls
        fold identically."""
        path = table(
            tmp_path,
            "fold.csv",
            [row().replace("\t", ",")],
            header=HEADER.replace("\t", ","),
        )
        parsed = read_typing_input([path])
        assert parsed.files[0].calls["marker_a"].identity_pct == 100.0

    def test_header_only_table_is_typed_no_data(self, tmp_path: Path) -> None:
        """Given a table whose screen found nothing (header only), When
        typed, Then exit 5 with TYPING_NO_DATA — the database to evaluate
        is only knowable from data rows."""
        path = table(tmp_path, "empty.tsv", [])
        result = runner.invoke(app, ["typing", str(path), "--datadir", str(typed_db)])
        assert result.exit_code == 5
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "TYPING_NO_DATA"

    def test_mixed_databases_are_a_mismatch(self, tmp_path: Path) -> None:
        """Given rows screening different databases, When typed, Then exit
        5 with DATABASE_MISMATCH listing the values."""
        path = table(
            tmp_path,
            "mixed.tsv",
            [row(database="markers"), row(file="t.fa", database="other")],
        )
        result = runner.invoke(app, ["typing", str(path), "--datadir", str(typed_db)])
        assert result.exit_code == 5
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "DATABASE_MISMATCH"
        assert envelope["context"]["databases"] == "markers, other"

    def test_missing_input_file(self, tmp_path: Path) -> None:
        """Given a path that does not exist, When typed, Then exit 5 with
        INPUT_NOT_FOUND."""
        result = runner.invoke(
            app, ["typing", str(tmp_path / "nope.tsv"), "--datadir", str(typed_db)]
        )
        assert result.exit_code == 5
        assert json.loads(result.stderr.splitlines()[-1])["code"] == "INPUT_NOT_FOUND"

    @pytest.mark.parametrize(
        ("rows", "detail"),
        [
            ([row()[:20]], "expected more columns than present"),
            ([row(identity="high")], "non-numeric %IDENTITY/%COVERAGE"),
        ],
    )
    def test_malformed_rows_fail_typed(self, tmp_path: Path, rows: list[str], detail: str) -> None:
        """Given a malformed data row, When typed, Then exit 5 with
        SCREEN_TABLE_MALFORMED naming the file and line."""
        path = table(tmp_path, "bad.tsv", rows)
        result = runner.invoke(app, ["typing", str(path), "--datadir", str(typed_db)])
        assert result.exit_code == 5
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "SCREEN_TABLE_MALFORMED"
        assert detail in envelope["message"]

    def test_missing_required_columns(self, tmp_path: Path) -> None:
        """Given a header without the DATABASE column, When typed, Then
        exit 5 with SCREEN_TABLE_MALFORMED naming the missing column."""
        path = table(
            tmp_path,
            "noheadercol.tsv",
            [row()],
            header=HEADER.replace("\tDATABASE", ""),
        )
        result = runner.invoke(app, ["typing", str(path), "--datadir", str(typed_db)])
        assert result.exit_code == 5
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert "DATABASE" in envelope["message"]


class TestCommand:
    def test_bare_invocation_at_a_terminal_prints_help(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Given no arguments at all with a terminal on stdin, When typing
        runs, Then the full help prints and exits 2 with no error envelope
        (the help-on-bare policy survives the pipe support; the tty probe
        is pinned because a CliRunner stdin is never a tty)."""
        monkeypatch.setattr(cmd_typing, "_stdin_is_tty", lambda: True)
        result = runner.invoke(app, ["typing"])
        assert result.exit_code == 2
        assert "Usage" in result.stdout
        assert "gapit.error" not in result.stderr

    def test_flags_without_input_reads_stdin(self, tmp_path: Path) -> None:
        """Given flags but no result file (an empty pipe on stdin), When
        typing runs, Then the typed no-data envelope fires — an omitted
        input now means stdin, not a usage error."""
        result = runner.invoke(app, ["typing", "--datadir", str(tmp_path)], input="")
        assert result.exit_code == 5
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "TYPING_NO_DATA"

    def test_csv_output_is_rejected(self, tmp_path: Path) -> None:
        """Given --format csv, When typing runs, Then a usage error names
        the available formats (tsv stays the joined default; csv input is
        still auto-detected)."""
        path = table(tmp_path, "ok.tsv", [row()])
        result = runner.invoke(app, ["typing", str(path), "--datadir", str(typed_db), "-f", "csv"])
        assert result.exit_code == 2
        assert "csv is not available" in json.loads(result.stderr.splitlines()[-1])["message"]

    def test_untyped_database_carries_no_scheme(self, tmp_path: Path, untyped_db: Path) -> None:
        """Given a table naming a database without typing.json, When typed,
        Then exit 4 with TYPING_NO_SCHEME naming the database."""
        path = table(tmp_path, "plain.tsv", [row()])
        result = runner.invoke(app, ["typing", str(path), "--datadir", str(untyped_db)])
        assert result.exit_code == 4
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "TYPING_NO_SCHEME"
        assert DB in envelope["message"]

    def test_cluster_database_typing_is_integrated(
        self, tmp_path: Path, tmp_path_factory: pytest.TempPathFactory
    ) -> None:
        """Given a table whose DATABASE names a cluster database, When
        typed, Then exit 4 with the typed error explaining cluster typing
        rides `gapit screen` itself."""
        datadir = tmp_path_factory.mktemp("typing_cmd_cluster")
        built = runner.invoke(
            app,
            [
                "db",
                "build",
                "loci",
                str(Path(__file__).parent / "data" / "cluster" / "screening.gbk"),
                "--datadir",
                str(datadir),
            ],
        )
        assert built.exit_code == 0, built.stderr
        path = table(tmp_path, "loci.tsv", [row(database="loci")])
        result = runner.invoke(app, ["typing", str(path), "--datadir", str(datadir)])
        assert result.exit_code == 4
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "TYPING_CLUSTER_DB"
        assert "gapit screen --db loci" in envelope["message"]

    def test_unknown_database_is_not_found(self, tmp_path: Path) -> None:
        """Given a table naming a database the datadir lacks, When typed,
        Then exit 4 with DATABASE_NOT_FOUND."""
        path = table(tmp_path, "ghost.tsv", [row(database="ghost")])
        result = runner.invoke(app, ["typing", str(path), "--datadir", str(tmp_path)])
        assert result.exit_code == 4
        assert json.loads(result.stderr.splitlines()[-1])["code"] == "DATABASE_NOT_FOUND"

    def test_output_writes_the_file_and_empties_stdout(
        self, tmp_path: Path, typed_db: Path
    ) -> None:
        """Given --output PATH, When typing runs, Then the results land in
        the file (header + rows), stdout carries no data, and the stderr
        note still prints."""
        path = table(tmp_path, "ok.tsv", [row()])
        target = tmp_path / "out.tsv"
        result = runner.invoke(
            app,
            [
                "typing",
                str(path),
                "--datadir",
                str(typed_db),
                "-o",
                str(target),
            ],
        )
        assert result.exit_code == 0, result.stderr
        assert result.stdout == ""
        assert "Typing 1 file(s)" in result.stderr
        written = target.read_text(encoding="utf-8")
        assert (
            written.splitlines()[0]
            == "FILE\tSCHEME\tPHENOTYPE\tCONFIDENCE\tSCORE\tRUNNER_UP\tNOTES"
        )

    def test_quiet_silences_the_stderr_note(self, tmp_path: Path, typed_db: Path) -> None:
        """Given --quiet, When typing runs, Then stderr is empty and the
        data still renders."""
        path = table(tmp_path, "ok.tsv", [row()])
        result = runner.invoke(app, ["typing", str(path), "--datadir", str(typed_db), "--quiet"])
        assert result.exit_code == 0
        assert result.stderr == ""
        assert result.stdout.startswith("FILE\tSCHEME\t")

    def test_short_and_long_forms_agree(self, tmp_path: Path, typed_db: Path) -> None:
        """Given the same invocation through the short and the long option
        spellings, When compared, Then the outputs are byte-identical (the
        single-dash alias table)."""
        path = table(tmp_path, "ok.tsv", [row()])
        short = runner.invoke(app, ["typing", str(path), "-D", str(typed_db), "-f", "json", "-q"])
        long = runner.invoke(
            app,
            [
                "typing",
                str(path),
                "--datadir",
                str(typed_db),
                "--format",
                "json",
                "--quiet",
            ],
        )
        assert short.exit_code == 0, short.stderr
        assert long.exit_code == 0, long.stderr
        assert short.stdout == long.stdout

    def test_json_document_shape(self, tmp_path: Path, typed_db: Path) -> None:
        """Given one screened table, When typed as JSON, Then the document
        is gapit.typing_result/1 with source paths, the db name, and one
        phenotypes object per FILE."""
        path = table(tmp_path, "shape.tsv", [row(), row(file="t.fa", gene="marker_b")])
        result = runner.invoke(
            app,
            ["typing", str(path), "--datadir", str(typed_db), "--format", "json", "--quiet"],
        )
        assert result.exit_code == 0, result.stderr
        document = json.loads(result.stdout)
        assert document["schema"] == "gapit.typing_result/1"
        assert document["tool"]["name"] == "gapit"
        assert document["source"] == [str(path)]
        assert document["db"] == DB
        assert [entry["file"] for entry in document["files"]] == ["s.fa", "t.fa"]
        assert set(document["files"][0]["phenotypes"]) == {"pathotype", "toxin"}

    def test_schema_registry_names_typing_result(self) -> None:
        """Given the schema registry, When the typing_result document is
        introspected, Then its JSON Schema names gapit.typing_result/1."""
        result = runner.invoke(app, ["schema", "typing_result"])
        assert result.exit_code == 0
        assert "gapit.typing_result/1" in result.stdout


class TestStdinInput:
    """The piped table: `gapit screen 1.fna --db NAME | gapit typing` types
    stdin with no file arguments; `-` is the explicit marker (it reads
    stdin even under a terminal, until EOF)."""

    def test_piped_table_types(self, typed_db: Path, tmp_path: Path) -> None:
        """Given a screen table on stdin with no file arguments, When
        typing runs, Then the TSV is byte-identical to the file-based run
        and the stderr note keeps its FILE-column count."""
        path = table(tmp_path, "ok.tsv", [row(), row(file="t.fa", gene="marker_b")])
        via_file = runner.invoke(app, ["typing", str(path), "--datadir", str(typed_db), "-q"])
        piped = runner.invoke(
            app,
            ["typing", "--datadir", str(typed_db), "-q"],
            input=path.read_text(encoding="utf-8"),
        )
        assert piped.exit_code == 0, piped.stderr
        assert piped.stdout == via_file.stdout
        note = runner.invoke(
            app,
            ["typing", "--datadir", str(typed_db)],
            input=path.read_text(encoding="utf-8"),
        )
        assert "Typing 2 file(s) against markers" in note.stderr

    def test_explicit_dash_matches_the_file_run(self, typed_db: Path, tmp_path: Path) -> None:
        """Given the `-` marker, When typing runs, Then stdin is the table
        and the output matches the file-based run byte for byte."""
        path = table(tmp_path, "ok.tsv", [row()])
        expected = runner.invoke(app, ["typing", str(path), "--datadir", str(typed_db), "-q"])
        dash = runner.invoke(
            app,
            ["typing", "-", "--datadir", str(typed_db), "-q"],
            input=path.read_text(encoding="utf-8"),
        )
        assert dash.exit_code == 0, dash.stderr
        assert dash.stdout == expected.stdout

    def test_dash_with_a_terminal_still_reads(
        self, typed_db: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Given `-` while stdin is a terminal, When typing runs, Then the
        marker still reads stdin until EOF — the tty guard only covers the
        argument-less invocation."""
        monkeypatch.setattr(cmd_typing, "_stdin_is_tty", lambda: True)
        result = runner.invoke(
            app,
            ["typing", "-", "--datadir", str(typed_db), "-q"],
            input="\n".join((HEADER, row())) + "\n",
        )
        assert result.exit_code == 0, result.stderr
        assert result.stdout.startswith("FILE\tSCHEME\t")

    def test_mixed_dash_and_files_is_a_usage_error(self, typed_db: Path, tmp_path: Path) -> None:
        """Given `-` alongside a file argument, When typing runs, Then a
        usage error demands one or the other (v1: either stdin or files,
        never both)."""
        path = table(tmp_path, "ok.tsv", [row()])
        result = runner.invoke(
            app,
            ["typing", "-", str(path), "--datadir", str(typed_db), "-q"],
            input="",
        )
        assert result.exit_code == 2
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "USAGE_ERROR"
        assert "not both" in envelope["message"]

    def test_empty_stdin_is_no_data(self, typed_db: Path) -> None:
        """Given an empty pipe, When typing runs, Then the existing
        TYPING_NO_DATA envelope fires, naming `-` as the source."""
        result = runner.invoke(app, ["typing", "--datadir", str(typed_db), "-q"], input="")
        assert result.exit_code == 5
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "TYPING_NO_DATA"
        assert envelope["message"].endswith("-")

    def test_json_source_names_stdin(self, typed_db: Path) -> None:
        """Given a piped table typed as JSON, When rendered, Then the
        source array names stdin as `-` (file runs keep their paths)."""
        result = runner.invoke(
            app,
            ["typing", "--datadir", str(typed_db), "-f", "json", "-q"],
            input="\n".join((HEADER, row())) + "\n",
        )
        assert result.exit_code == 0, result.stderr
        assert json.loads(result.stdout)["source"] == ["-"]

    def test_malformed_piped_row_names_stdin(self, typed_db: Path) -> None:
        """Given a malformed piped row, When typing runs, Then the
        SCREEN_TABLE_MALFORMED envelope names `-` and the line number."""
        result = runner.invoke(
            app,
            ["typing", "--datadir", str(typed_db), "-q"],
            input="\n".join((HEADER, row()[:20])) + "\n",
        )
        assert result.exit_code == 5
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "SCREEN_TABLE_MALFORMED"
        assert "(-:2)" in envelope["message"]

    def test_file_dash_and_pipe_produce_identical_tsv(self, typed_db: Path, tmp_path: Path) -> None:
        """Given one real screen run's table delivered three ways — written
        with `-o` then typed by path, typed through the `-` marker, and
        piped bare (`screen | typing`) — When compared, Then all three
        typing outputs are byte-identical TSV."""
        samples = [str(DATA / "exact.fa"), str(DATA / "partial.fa"), str(DATA / "mutated.fa")]
        table_path = tmp_path / "screen.tsv"
        to_file = runner.invoke(
            app,
            [
                "screen",
                *samples,
                "--db",
                DB,
                "--datadir",
                str(typed_db),
                "-o",
                str(table_path),
                "-q",
            ],
        )
        assert to_file.exit_code == 0, to_file.stderr
        on_stdout = runner.invoke(
            app, ["screen", *samples, "--db", DB, "--datadir", str(typed_db), "-q"]
        )
        assert on_stdout.exit_code == 0, on_stdout.stderr

        by_path = runner.invoke(app, ["typing", str(table_path), "--datadir", str(typed_db), "-q"])
        by_dash = runner.invoke(
            app,
            ["typing", "-", "--datadir", str(typed_db), "-q"],
            input=table_path.read_text(encoding="utf-8"),
        )
        by_pipe = runner.invoke(
            app, ["typing", "--datadir", str(typed_db), "-q"], input=on_stdout.stdout
        )
        assert by_path.exit_code == 0, by_path.stderr
        assert by_dash.exit_code == 0, by_dash.stderr
        assert by_pipe.exit_code == 0, by_pipe.stderr
        assert by_path.stdout == by_dash.stdout == by_pipe.stdout
        assert by_path.stdout.startswith("FILE\tSCHEME\t")


@pytest.fixture(scope="module")
def markers_document(typed_db: Path):
    database = next(entry for entry in discover_databases(typed_db) if entry.name == DB)
    typing_document = load_typing(database)
    assert typing_document is not None
    return typing_document


class TestEquivalence:
    """The two-stage pipeline against the inline engine: for each fixture
    sample, `gapit typing` over the screened table must produce exactly the
    calls evaluate_gene_calls produces over the screened Report's folded
    best hits (the pre-purification inline semantics)."""

    SAMPLES = ("exact", "partial", "mutated")

    def _inline_calls(self, hits: list[dict[str, Any]], document: Any) -> dict[str, Any]:
        best: dict[str, dict[str, Any]] = {}
        for hit in hits:
            current = best.get(hit["gene"])
            if current is None or (hit["identity_pct"], hit["coverage_pct"]) > (
                current["identity_pct"],
                current["coverage_pct"],
            ):
                best[hit["gene"]] = hit
        calls = {
            gene: GeneCall(
                gene_id=gene,
                start=0,
                end=0,
                strand="+",
                coverage_pct=hit["coverage_pct"],
                identity_pct=hit["identity_pct"],
                verdict="present",
            )
            for gene, hit in best.items()
        }
        return {
            scheme: call.model_dump(mode="json")
            for scheme, call in evaluate_gene_calls(calls, document).items()
        }

    def test_typing_output_equals_the_inline_engine(
        self, typed_db: Path, markers_document: Any, tmp_path: Path
    ) -> None:
        """Given the three marker samples, When screened as JSON (for the
        engine) and as a TSV table (for the command), Then every FILE's
        typing calls are identical to the inline engine's calls."""
        screen_json = runner.invoke(
            app,
            [
                "screen",
                *(str(DATA / f"{sample}.fa") for sample in self.SAMPLES),
                "--db",
                DB,
                "--datadir",
                str(typed_db),
                "--format",
                "json",
                "--quiet",
            ],
        )
        assert screen_json.exit_code == 0, screen_json.stderr
        table_path = tmp_path / "screen.tsv"
        screened = runner.invoke(
            app,
            [
                "screen",
                *(str(DATA / f"{sample}.fa") for sample in self.SAMPLES),
                "--db",
                DB,
                "--datadir",
                str(typed_db),
                "-o",
                str(table_path),
                "--quiet",
            ],
        )
        assert screened.exit_code == 0, screened.stderr
        result = runner.invoke(
            app,
            ["typing", str(table_path), "--datadir", str(typed_db), "--format", "json", "--quiet"],
        )
        assert result.exit_code == 0, result.stderr
        typed = {entry["file"]: entry["phenotypes"] for entry in json.loads(result.stdout)["files"]}
        for entry in json.loads(screen_json.stdout)["files"]:
            assert typed[entry["file"]] == self._inline_calls(entry["hits"], markers_document)


class TestGoldens:
    def test_json_and_md_goldens_at_pinned_now(
        self,
        typed_db: Path,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Given the markers samples typed through the direct renderers at
        the pinned now, Then each output is byte-identical to its committed
        golden (source pinned to the literal table name)."""
        monkeypatch.chdir(DATA)
        table_path = tmp_path / "screen.tsv"
        screened = runner.invoke(
            app,
            [
                "screen",
                "exact.fa",
                "partial.fa",
                "mutated.fa",
                "--db",
                DB,
                "--datadir",
                str(typed_db),
                "-o",
                str(table_path),
                "--quiet",
            ],
        )
        assert screened.exit_code == 0, screened.stderr
        database = next(entry for entry in discover_databases(typed_db) if entry.name == DB)
        document = load_typing(database)
        assert document is not None
        results = [
            TypingFileResult(file=entry.file, phenotypes=evaluate_gene_calls(entry.calls, document))
            for entry in read_typing_input([table_path]).files
        ]
        assert render_typing_result_json(["screen.tsv"], DB, results, now=PINNED_NOW) == (
            GOLDEN / "typing_markers.json"
        ).read_text(encoding="utf-8")
        assert render_typing_md(["screen.tsv"], DB, results, now=PINNED_NOW) == (
            GOLDEN / "typing_markers.md"
        ).read_text(encoding="utf-8")
