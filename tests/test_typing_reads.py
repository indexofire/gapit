"""Reads-mode input for `gapit typing`: #SAMPLE/BREADTH% tables.

The second stage originally read only contig screen tables. A reads screen
(`gapit screen --r1/--r2`) renders the same gene calls under different
column names, so typing must fold those rows too: #SAMPLE is the FILE key,
BREADTH% is the coverage metric, PRESENT=no rows are skipped, and reads/2's
opt-in IDENTITY% supplies the identity when present.
"""

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from gapit.cli import app
from gapit.fasta import iter_fasta
from gapit.typing_input import read_typing_input

DATA = Path(__file__).parent / "data" / "typing"
DB = "markers"
READS_HEADER = "#SAMPLE\tGENE\tBREADTH%\tDEPTH\tREADS\tPRESENT\tDATABASE\tACCESSION\tPRODUCT"
READS2_HEADER = f"{READS_HEADER}\tIDENTITY%"

runner = CliRunner()


def reads_table(tmp_path: Path, rows: list[str], *, header: str = READS_HEADER) -> Path:
    """A reads result table from data lines (tab-separated)."""
    path = tmp_path / "reads.tsv"
    path.write_text("\n".join((header, *rows)) + "\n", encoding="utf-8")
    return path


def reads_row(
    *,
    sample: str = "reads.fq",
    gene: str = "marker_a",
    breadth: str = "99.50",
    present: str = "yes",
    database: str = DB,
    identity: str | None = None,
) -> str:
    cells = [sample, gene, breadth, "12.00", "8", present, database, "acc", "prod"]
    if identity is not None:
        cells.append(identity)
    return "\t".join(cells)


@pytest.fixture(scope="module")
def typed_db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    datadir = tmp_path_factory.mktemp("typing_reads")
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


class TestParsing:
    def test_reads1_breadth_is_coverage_and_identity_is_unknown(self, tmp_path: Path) -> None:
        """Given a reads/1 row, When parsed, Then #SAMPLE is the FILE key,
        BREADTH% is the coverage, and the missing identity takes the
        reads/1 sentinel (the format cannot gate identity)."""
        path = reads_table(tmp_path, [reads_row(breadth="97.25")])
        (entry,) = read_typing_input([path]).files
        assert entry.file == "reads.fq"
        call = entry.calls["marker_a"]
        assert call.coverage_pct == 97.25
        assert call.identity_pct == 100.0
        assert call.verdict == "present"

    def test_reads2_identity_column_is_used(self, tmp_path: Path) -> None:
        """Given a reads/2 row, When parsed, Then IDENTITY% supplies the
        identity instead of the reads/1 sentinel."""
        path = reads_table(
            tmp_path,
            [reads_row(identity="98.75")],
            header=READS2_HEADER,
        )
        (entry,) = read_typing_input([path]).files
        call = entry.calls["marker_a"]
        assert call.identity_pct == 98.75
        assert call.coverage_pct == 99.5

    def test_absent_rows_are_skipped(self, tmp_path: Path) -> None:
        """Given --all-genes output, When parsed, Then PRESENT=no rows never
        become gene calls (typing scores presence, not absence rows)."""
        path = reads_table(
            tmp_path,
            [
                reads_row(gene="marker_a", present="yes"),
                reads_row(gene="marker_b", present="no"),
            ],
        )
        (entry,) = read_typing_input([path]).files
        assert set(entry.calls) == {"marker_a"}

    def test_missing_breadth_names_the_reads_column(self, tmp_path: Path, typed_db: Path) -> None:
        """Given a reads-shaped header without BREADTH%, When typed, Then
        the malformed envelope names BREADTH%, not the contig columns."""
        header = "#SAMPLE\tGENE\tDEPTH\tREADS\tPRESENT\tDATABASE\tACCESSION\tPRODUCT"
        path = reads_table(tmp_path, [], header=header)
        result = runner.invoke(app, ["typing", str(path), "--datadir", str(typed_db), "-q"])
        assert result.exit_code == 5
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "SCREEN_TABLE_MALFORMED"
        assert "BREADTH%" in envelope["message"]


class TestCommand:
    def test_reads_table_types(self, tmp_path: Path, typed_db: Path) -> None:
        """Given a reads result table, When typed, Then the sample gets its
        phenotype call from breadth-based presence."""
        path = reads_table(tmp_path, [reads_row()])
        result = runner.invoke(app, ["typing", str(path), "--datadir", str(typed_db), "-q"])
        assert result.exit_code == 0, result.stderr
        lines = result.stdout.splitlines()
        assert lines[0].startswith("FILE\tSCHEME\t")
        assert "reads.fq" in lines[1]
        assert "marker_a" in lines[1]
        assert "EPEC" in lines[1]

    def test_fastq_screen_pipes_into_typing(self, tmp_path: Path, typed_db: Path) -> None:
        """Given real FASTQ screened through minimap2, When its table is
        piped into typing, Then the end-to-end pipeline designates the
        sample (the issue's exact workflow)."""
        reads = tmp_path / "marker_a.fq"
        record = next(
            record for record in iter_fasta(DATA / "markers.fa") if record.id == "marker_a"
        )
        with reads.open("w", encoding="utf-8") as handle:
            for index, start in enumerate(range(0, len(record.sequence), 150), start=1):
                fragment = record.sequence[start : start + 150]
                if len(fragment) < 50:
                    continue
                handle.write(f"@read{index}\n{fragment}\n+\n{'I' * len(fragment)}\n")

        screened = runner.invoke(
            app,
            [
                "screen",
                "--r1",
                str(reads),
                "--db",
                DB,
                "--datadir",
                str(typed_db),
                "--quiet",
            ],
        )
        assert screened.exit_code == 0, screened.stderr
        assert screened.stdout.startswith("#SAMPLE\t")

        typed = runner.invoke(
            app,
            ["typing", "--datadir", str(typed_db), "--quiet"],
            input=screened.stdout,
        )
        assert typed.exit_code == 0, typed.stderr
        assert "marker_a" in typed.stdout
        assert "EPEC" in typed.stdout
