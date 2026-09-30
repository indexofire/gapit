"""Offline self-test of scripts/cluster_calibration.py (developer tool).

Loads the script by file path (scripts/ is not a package) behind a
runtime-checked Protocol, builds the typed fixture db, writes a labels TSV,
and runs the harness end to end against the real engine — asserting the
three report blocks on known calls (the integration seam the future
training loop will ride).
"""

import importlib.util
from pathlib import Path
from typing import Protocol, runtime_checkable

import pytest
from typer.testing import CliRunner

from gapit.cli import app
from gapit.fasta import iter_fasta

SCRIPT = Path(__file__).parent.parent / "scripts" / "cluster_calibration.py"
DATA = Path(__file__).parent / "data" / "cluster"


@runtime_checkable
class Harness(Protocol):
    """The harness surface these tests drive (main + the labels reader)."""

    def main(self, argv: list[str] | None = ...) -> int: ...

    def read_labels(self, path: Path) -> list[tuple[str, str]]: ...


def _load_harness() -> Harness:
    spec = importlib.util.spec_from_file_location("cluster_calibration", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert isinstance(module, Harness)
    return module


@pytest.fixture()
def harness() -> Harness:
    return _load_harness()


@pytest.fixture()
def typed_env(tmp_path: Path) -> tuple[Path, dict[str, Path]]:
    """The typed fixture db plus exact/both/unrelated samples."""
    runner = CliRunner()
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    result = runner.invoke(
        app,
        [
            "db",
            "build",
            "cps",
            str(DATA / "screening.gbk"),
            "--datadir",
            str(datadir),
            "--typing",
            str(DATA / "typing_screen.json"),
        ],
    )
    assert result.exit_code == 0, result.stderr
    loci = {record.id: record.sequence for record in iter_fasta(datadir / "cps" / "sequences")}
    samples = tmp_path / "samples"
    samples.mkdir()
    (samples / "exact.fa").write_text(f">ctg_a\n{loci['locusA']}\n", encoding="utf-8")
    (samples / "both.fa").write_text(
        f">ctg_a\n{loci['locusA']}\n>ctg_b\n{loci['locusB']}\n", encoding="utf-8"
    )
    (samples / "unrelated.fa").write_text(">ctg\n" + "ACGT" * 250 + "\n", encoding="utf-8")
    return datadir, {path.stem: path for path in sorted(samples.glob("*.fa"))}


def test_calibration_reports_distribution_matrix_and_divergences(
    harness: Harness,
    typed_env: tuple[Path, dict[str, Path]],
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Given the typed fixture db and labels (exact=K101, both=K101,
    unrelated=K101), When the harness runs, Then the distribution block
    carries K101's score stats, the matrix counts the one agreeing call,
    and the divergences list the ambiguous and no-locus samples."""
    datadir, samples = typed_env
    labels = tmp_path / "labels.tsv"
    labels.write_text(
        "\n".join(
            [
                "# sample\texpected",
                f"{samples['exact']}\tK101",
                f"{samples['both']}\tK101",
                f"{samples['unrelated']}\tK101",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    exit_code = harness.main(["cps", "--datadir", str(datadir), str(labels)])

    assert exit_code == 0
    out = capsys.readouterr().out
    assert "expected\tn\tmin\tmedian\tmax" in out
    # the K101 rule scores 1.0 (exact), 1.0 (both), 0.0 (unrelated: require_any unsatisfied)
    assert "K101\t3\t0.0000\t1.0000\t1.0000" in out
    assert "agreement matrix" in out
    matrix_rows = [line for line in out.splitlines() if line.startswith("K101\t")]
    assert matrix_rows and matrix_rows[0].count("\t1") >= 1  # one K101->K101 cell
    assert "divergences (2 of 3)" in out
    assert f"{samples['both']}\tK101\t-\t1.0000" in out
    assert f"{samples['unrelated']}\tK101\t(no call)\t0.0000" in out


def test_calibration_refuses_untyped_database(
    harness: Harness, tmp_path: Path, typed_env: tuple[Path, dict[str, Path]]
) -> None:
    """Given the same fixture db rebuilt WITHOUT a typing spec, When the
    harness runs, Then it refuses with the typed-database message."""
    datadir, samples = typed_env
    (datadir / "cps" / "typing.json").unlink()
    labels = tmp_path / "labels.tsv"
    labels.write_text(f"{samples['exact']}\tK101\n", encoding="utf-8")
    with pytest.raises(SystemExit) as raised:
        harness.main(["cps", "--datadir", str(datadir), str(labels)])
    assert "no typing.json" in str(raised.value)


def test_calibration_json_roundtrip_of_labels(
    harness: Harness, typed_env: tuple[Path, dict[str, Path]], tmp_path: Path
) -> None:
    """Given labels with comments, blanks, and a malformed row, When read,
    Then only well-formed rows survive (the read_labels contract)."""
    _, samples = typed_env
    labels = tmp_path / "labels.tsv"
    labels.write_text(
        "\n".join(
            [
                "# comment",
                "",
                f"{samples['exact']}\tK101",
                "no-tab-line",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    rows = harness.read_labels(labels)
    assert rows == [(str(samples["exact"]), "K101")]
