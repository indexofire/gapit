"""Bundled Shigella marker panel (first serotype-DB slice)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TextIO

import pytest
from click.testing import CliRunner

from gspec.bundled import find_bundled
from gspec.cli import app
from gspec.dbcodec import decode_seqid
from gspec.fasta import iter_fasta

runner = CliRunner()


def _marker_sequences(marker: str) -> list[str]:
    bundled = find_bundled("shi_sero")
    assert bundled is not None
    return [
        record.sequence
        for record in iter_fasta(bundled.sequences_path)
        if decode_seqid(record.id, "shi_sero").gene == marker
    ]


def _marker_sequence(marker: str) -> str:
    return _marker_sequences(marker)[0]


def test_bundled_shi_sero_has_marker_metadata() -> None:
    """Given the bundled panel, When inspected, Then it carries the 95
    ShigaTyper markers (ipaH_c split into two control fragments) and the
    O-antigen/POAC/virulence roles the profile resolver consumes."""
    bundled = find_bundled("shi_sero")
    assert bundled is not None
    markers = json.loads((bundled.path / "markers.json").read_text(encoding="utf-8"))
    assert markers["schema"] == "gspec.shigella_markers/1"
    assert len(markers["markers"]) == 96
    categories = {entry["category"] for entry in markers["markers"]}
    assert {"o_antigen", "poac", "virulence"} <= categories
    groups = {entry["serotype_group"] for entry in markers["markers"] if entry["serotype_group"]}
    assert {"SD1", "SB2", "SF1-5", "SF6", "SS"} <= groups
    controls = [entry for entry in markers["markers"] if entry["gene"] == "ipaH_c"]
    assert {entry["length"] for entry in controls} == {640, 140}
    assert {entry["derived_from"] for entry in controls} == {"ipaH_c:1-640", "ipaH_c:641-780"}


def test_bundled_shi_sero_marker_screens(tmp_path: Path) -> None:
    """Given a deposited marker, When screened against the bundled panel,
    Then the marker is detected through the standard gene engine."""
    sample = tmp_path / "sd1.fna"
    sample.write_text(f">sd1\n{_marker_sequence('Sd1_wzx')}\n", encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    result = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            "shi_sero",
            "--datadir",
            str(datadir),
            "--quiet",
        ],
    )
    assert result.exit_code == 0, result.stderr
    assert "Sd1_wzx" in result.stdout
    assert "shi_sero" in result.stdout
    assert (datadir / "shi_sero" / "markers.json").is_file()
    assert (datadir / "shi_sero" / "serotype.json").is_file()


def test_bundled_shi_sero_has_serotype_profile() -> None:
    """Given the bundle, When the profile is loaded, Then it is the first
    gspec.serotype/1 rule set with ipaH as the control, and the SB6/SB10
    rules carry the documented wbaM 5'-junction gate."""
    bundled = find_bundled("shi_sero")
    assert bundled is not None
    document = json.loads((bundled.path / "serotype.json").read_text(encoding="utf-8"))
    assert document["schema"] == "gspec.serotype/1"
    assert document["control_gene"] == "ipaH_c"
    assert len(document["rules"]) == 63
    rules = {rule["phenotype"]: rule for rule in document["rules"]}
    assert rules["S. boydii 10"]["max_alignment_start_pct"] == {"wbaM": 20.0}
    assert rules["S. boydii 6"]["min_alignment_start_pct"] == {"wbaM": 20.0}
    assert rules["S. boydii 6/10"]["excludes"] == ["wbaM"]
    screen = json.loads((bundled.path / "screen.json").read_text(encoding="utf-8"))
    assert screen["contigs"]["merge_fragments"] is True


def test_bundled_shi_sero_types_sd1(tmp_path: Path) -> None:
    """Given ipaH + Sd1 wzx/rfp markers, When typed, Then the serotype
    profile resolves S. dysenteriae 1 through the standard typing table."""

    def _write_marker(handle: TextIO, marker: str) -> None:
        handle.write(f">{marker}\n{_marker_sequence(marker)}\n")

    sample = tmp_path / "sd1.fna"
    with sample.open("w", encoding="utf-8") as handle:
        for marker in ("ipaH_c", "Sd1_wzx", "Sd1_rfp"):
            _write_marker(handle, marker)
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    screened = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            "shi_sero",
            "--datadir",
            str(datadir),
            "--quiet",
            "-o",
            str(tmp_path / "screen.tsv"),
        ],
    )
    assert screened.exit_code == 0, screened.stderr
    typed = runner.invoke(
        app,
        [
            "typing",
            str(tmp_path / "screen.tsv"),
            "--datadir",
            str(datadir),
            "--quiet",
        ],
    )
    assert typed.exit_code == 0, typed.stderr
    rows = typed.stdout.splitlines()
    assert rows[0].startswith("FILE\tSCHEME\tPHENOTYPE")
    assert rows[1].split("\t")[:3] == [str(sample), "serotype", "S. dysenteriae 1"]


@pytest.mark.parametrize(
    ("start", "end", "expected"),
    [
        (0, 831, "S. boydii 10"),
        (252, None, "S. boydii 6"),
    ],
)
def test_bundled_shi_sero_separates_sb6_sb10_by_wbam_start(
    tmp_path: Path, start: int, end: int | None, expected: str
) -> None:
    """Given a 3'-truncated wbaM, When typed, Then the 5'-intact alignment
    calls SB10; given a 5'-truncated wbaM (the IS629 junction), SB6."""
    sample = tmp_path / "sb.fna"
    sample.write_text(
        "".join(
            f">{marker}\n{sequence}\n"
            for marker, sequence in (
                ("ipaH_c", _marker_sequence("ipaH_c")),
                ("Sb6_wzx", _marker_sequence("Sb6_wzx")),
                ("wbaM", _marker_sequence("wbaM")[start:end]),
            )
        ),
        encoding="utf-8",
    )
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    screened = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            "shi_sero",
            "--datadir",
            str(datadir),
            "--quiet",
            "-o",
            str(tmp_path / "screen.tsv"),
        ],
    )
    assert screened.exit_code == 0, screened.stderr
    typed = runner.invoke(
        app,
        ["typing", str(tmp_path / "screen.tsv"), "--datadir", str(datadir), "--quiet"],
    )
    assert typed.exit_code == 0, typed.stderr
    assert typed.stdout.splitlines()[1].split("\t")[2] == expected


@pytest.mark.parametrize("control_index", [0, 1])
def test_bundled_shi_sero_uses_split_ipaH_control_fragments(
    tmp_path: Path, control_index: int
) -> None:
    """Given either half of the split ipaH control, When typed, Then the
    serotype resolves through the standard control gate."""
    control = _marker_sequences("ipaH_c")[control_index]
    assert len(control) in (640, 140)
    sample = tmp_path / "sd1_short_control.fna"
    sample.write_text(
        "".join(
            f">{marker}\n{sequence}\n"
            for marker, sequence in (
                ("ipaH_c", control),
                ("Sd1_wzx", _marker_sequence("Sd1_wzx")),
                ("Sd1_rfp", _marker_sequence("Sd1_rfp")),
            )
        ),
        encoding="utf-8",
    )
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    screened = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            "shi_sero",
            "--datadir",
            str(datadir),
            "--quiet",
            "-o",
            str(tmp_path / "screen.tsv"),
        ],
    )
    assert screened.exit_code == 0, screened.stderr
    typed = runner.invoke(
        app,
        ["typing", str(tmp_path / "screen.tsv"), "--datadir", str(datadir), "--quiet"],
    )
    assert typed.exit_code == 0, typed.stderr
    row = typed.stdout.splitlines()[1].split("\t")
    assert row[2] == "S. dysenteriae 1"


def test_bundled_shi_sero_without_control_is_not_shigella(tmp_path: Path) -> None:
    """Given a marker but no ipaH, When typed, Then the control gate falls
    back instead of reporting a serotype."""
    sample = tmp_path / "marker.fna"
    sample.write_text(f">sd1\n{_marker_sequence('Sd1_wzx')}\n", encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    screened = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            "shi_sero",
            "--datadir",
            str(datadir),
            "--quiet",
            "-o",
            str(tmp_path / "screen.tsv"),
        ],
    )
    assert screened.exit_code == 0, screened.stderr
    typed = runner.invoke(
        app,
        ["typing", str(tmp_path / "screen.tsv"), "--datadir", str(datadir), "--quiet"],
    )
    assert typed.exit_code == 0, typed.stderr
    row = typed.stdout.splitlines()[1].split("\t")
    assert row[1:4] == ["serotype", "Not Shigella/EIEC", "Sd1_wzx"]


def test_db_info_exposes_the_serotype_schema(tmp_path: Path) -> None:
    """Given the bundled marker panel, When db info runs before install,
    Then it advertises gspec.serotype/1 and the typing workflow."""
    result = runner.invoke(
        app,
        [
            "db",
            "info",
            "-d",
            "shi_sero",
            "--datadir",
            str(tmp_path / "datadir"),
            "--json",
        ],
    )
    assert result.exit_code == 0, result.stderr
    document = json.loads(result.stdout)
    assert document["typing_schema"] == "gspec.serotype/1"
    assert any("gspec typing" in line for line in document["usage"])
