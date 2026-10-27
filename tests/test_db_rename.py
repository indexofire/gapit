"""`gapit db rename`: re-homing an installed database under a new name.

Gene databases carry the name in every record and header, so a rename must
rebuild them; cluster databases carry it only in the manifest. Both paths are
staged, so a failed rename leaves the source intact.
"""

import json
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from gapit.bundled import find_bundled, materialize_bundled
from gapit.cli import app
from gapit.dbcodec import decode_seqid
from gapit.fasta import iter_fasta
from gapit.records import read_manifest

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"
READS = Path(__file__).parent / "data" / "reads"

runner = CliRunner()


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    """The committed fixture datadir (tinyamr built)."""
    target = tmp_path / "datadir"
    shutil.copytree(FIXTURE_DB_DIR, target)
    return target


def _gene_db_from_bundle(name: str, datadir: Path, target: str) -> None:
    """Materialize a bundled gene database under ``target`` (real records)."""
    bundled = find_bundled(name)
    assert bundled is not None
    command = ["db", "build", target, str(bundled.sequences_path), "--datadir", str(datadir)]
    if bundled.typing_path is not None:
        command.extend(["--typing", str(bundled.typing_path)])
    built = runner.invoke(app, command)
    assert built.exit_code == 0, built.stderr


def test_rename_gene_database_rewrites_records_headers_and_manifest(
    datadir: Path,
) -> None:
    """Given an installed gene database, When renamed, Then the receipt names
    the new database, the old directory is gone, and every record, header, and
    manifest carries the new name while the sequences survive byte-for-byte."""
    _gene_db_from_bundle("lm_doumith", datadir, "old_lm")
    before = {
        (decode_seqid(record.id, default_db="old_lm").gene, record.sequence)
        for record in iter_fasta(datadir / "old_lm" / "sequences")
    }

    result = runner.invoke(app, ["db", "rename", "old_lm", "lm_panel", "--datadir", str(datadir)])

    assert result.exit_code == 0, result.stderr
    receipt = json.loads(result.stdout)
    assert receipt == {
        "db": "lm_panel",
        "renamed_from": "old_lm",
        "records": 5,
        "kind": "gene",
        "destination": str(datadir / "lm_panel"),
    }
    assert not (datadir / "old_lm").exists()
    manifest = read_manifest(datadir / "lm_panel" / "gapit-manifest.json")
    assert manifest.name == "lm_panel"
    assert manifest.n_records == 5
    records = [
        json.loads(line)
        for line in (datadir / "lm_panel" / "records.jsonl")
        .read_text(encoding="utf-8")
        .splitlines()
    ]
    assert {record["db"] for record in records} == {"lm_panel"}
    after = {
        (decode_seqid(record.id, default_db="lm_panel").gene, record.sequence)
        for record in iter_fasta(datadir / "lm_panel" / "sequences")
    }
    assert after == before
    for record in iter_fasta(datadir / "lm_panel" / "sequences"):
        assert decode_seqid(record.id, default_db="lm_panel").database == "lm_panel"
    # The typing sidecar travels along, so the renamed db still designates.
    assert (datadir / "lm_panel" / "typing.json").is_file()


def test_renamed_gene_database_still_screens_and_types(datadir: Path) -> None:
    """Given a renamed typed gene database, When the two-stage pipeline runs,
    Then the screen names the new database and typing still designates."""
    _gene_db_from_bundle("lm_doumith", datadir, "old_lm")
    renamed = runner.invoke(app, ["db", "rename", "old_lm", "lm_panel", "--datadir", str(datadir)])
    assert renamed.exit_code == 0, renamed.stderr
    bundled = find_bundled("lm_doumith")
    assert bundled is not None
    sample = datadir.parent / "sample.fa"
    sample.write_text(
        "".join(
            f">{decode_seqid(record.id, default_db='lm_panel').gene}\n{record.sequence}\n"
            for record in iter_fasta(bundled.sequences_path)
        ),
        encoding="utf-8",
    )
    table = datadir.parent / "sample.tsv"
    screened = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            "lm_panel",
            "--datadir",
            str(datadir),
            "--output",
            str(table),
            "--quiet",
        ],
    )
    assert screened.exit_code == 0, screened.stderr
    assert {line.split("\t")[11] for line in table.read_text().splitlines()[1:]} == {"lm_panel"}
    typed = runner.invoke(
        app,
        ["typing", str(table), "--datadir", str(datadir), "--format", "json", "--quiet"],
    )
    assert typed.exit_code == 0, typed.stderr
    assert json.loads(typed.stdout)["db"] == "lm_panel"


def test_rename_cluster_database_copies_artifacts_and_restamps_manifest(
    datadir: Path,
) -> None:
    """Given an installed cluster database, When renamed, Then the locus FASTA
    and feature table travel verbatim (their headers carry no db name), the
    manifest and its record count move, and screening still works."""
    materialize_bundled("vp_ok", datadir, quiet=True)
    sequences_bytes = (datadir / "vp_ok" / "sequences").read_bytes()

    result = runner.invoke(app, ["db", "rename", "vp_ok", "vpa_ok", "--datadir", str(datadir)])

    assert result.exit_code == 0, result.stderr
    receipt = json.loads(result.stdout)
    assert (receipt["db"], receipt["kind"], receipt["records"]) == ("vpa_ok", "cluster", 195)
    assert not (datadir / "vp_ok").exists()
    assert (datadir / "vpa_ok" / "sequences").read_bytes() == sequences_bytes
    assert (datadir / "vpa_ok" / "features.json").is_file()
    manifest = read_manifest(datadir / "vpa_ok" / "gapit-manifest.json")
    assert (manifest.name, manifest.kind, manifest.header_format) == ("vpa_ok", "cluster", "plain")
    assert (datadir / "vpa_ok" / "sequences.nin").is_file()


def test_rename_refuses_a_missing_source(datadir: Path) -> None:
    """Given a name that is not installed, Then exit 4 DATABASE_NOT_FOUND."""
    result = runner.invoke(app, ["db", "rename", "nope", "other", "--datadir", str(datadir)])
    assert result.exit_code == 4
    envelope = json.loads(result.stderr.splitlines()[-1])
    assert envelope["code"] == "DATABASE_NOT_FOUND"


def test_rename_refuses_an_existing_target_unless_forced(datadir: Path) -> None:
    """Given an installed target, Then exit 4 DB_ALREADY_EXISTS; --force
    replaces it and still leaves exactly one copy of each name."""
    _gene_db_from_bundle("lm_doumith", datadir, "one")
    _gene_db_from_bundle("lm_doumith", datadir, "two")

    refused = runner.invoke(app, ["db", "rename", "one", "two", "--datadir", str(datadir)])
    assert refused.exit_code == 4
    assert json.loads(refused.stderr.splitlines()[-1])["code"] == "DB_ALREADY_EXISTS"
    assert (datadir / "one").is_dir()

    forced = runner.invoke(
        app, ["db", "rename", "one", "two", "--datadir", str(datadir), "--force"]
    )
    assert forced.exit_code == 0, forced.stderr
    assert not (datadir / "one").exists()
    assert read_manifest(datadir / "two" / "gapit-manifest.json").name == "two"


def test_rename_rejects_bad_names_and_self_renames(datadir: Path) -> None:
    """Given an escaping target name or an identical pair, Then the call fails
    as a usage/DB error instead of touching the datadir."""
    escaping = runner.invoke(app, ["db", "rename", "tinyamr", "../x", "--datadir", str(datadir)])
    assert escaping.exit_code == 2
    same = runner.invoke(app, ["db", "rename", "tinyamr", "tinyamr", "--datadir", str(datadir)])
    assert same.exit_code == 4
    assert (datadir / "tinyamr" / "sequences").is_file()
