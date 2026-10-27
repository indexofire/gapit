"""`gapit db build NAME` (no INPUT): re-project a gene database from its records.

`records.jsonl` is the editable truth, so gene curation (renaming a gene,
merging alleles under one name) is a records edit plus a rebuild. The rebuild
must be staged (a broken edit leaves the database usable) and must refuse
cluster databases, whose artifacts are generated from a source GBK/GFF.
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

runner = CliRunner()


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    target = tmp_path / "datadir"
    shutil.copytree(FIXTURE_DB_DIR, target)
    return target


def _materialize(name: str, target: str, datadir: Path) -> None:
    bundled = find_bundled(name)
    assert bundled is not None
    built = runner.invoke(
        app, ["db", "build", target, str(bundled.sequences_path), "--datadir", str(datadir)]
    )
    assert built.exit_code == 0, built.stderr


def test_rebuild_republishes_an_edited_gene_name(datadir: Path) -> None:
    """Given a database whose records were edited to rename a gene, When
    `db build NAME` runs, Then the sequences headers, the manifest count, and the
    BLAST index all follow the edit."""
    _materialize("lm_doumith", "lm", datadir)
    records_path = datadir / "lm" / "records.jsonl"
    records = [json.loads(line) for line in records_path.read_text(encoding="utf-8").splitlines()]
    renamed = [
        {**record, "gene": "prs_renamed"} if record["gene"] == "prs" else record
        for record in records
    ]
    records_path.write_text(
        "".join(json.dumps(record) + "\n" for record in renamed), encoding="utf-8"
    )

    result = runner.invoke(app, ["db", "build", "lm", "--datadir", str(datadir)])

    assert result.exit_code == 0, result.stderr
    receipt = json.loads(result.stdout)
    assert receipt == {
        "db": "lm",
        "records": 5,
        "dbtype": "nucl",
        "destination": str(datadir / "lm"),
    }
    genes = {
        decode_seqid(record.id, default_db="lm").gene
        for record in iter_fasta(datadir / "lm" / "sequences")
    }
    assert "prs_renamed" in genes and "prs" not in genes
    assert read_manifest(datadir / "lm" / "gapit-manifest.json").n_records == 5
    assert (datadir / "lm" / "sequences.nin").is_file()


def test_rebuild_merges_two_alleles_under_one_gene_name(datadir: Path) -> None:
    """Given two records renamed to the same gene (alleles of one locus), When
    the database is rebuilt, Then both records survive (the duplicate dropper
    only drops identical gene+sequence pairs) and a screen reports the single
    shared name."""
    bundled = find_bundled("lm_doumith")
    assert bundled is not None
    fasta = datadir.parent / "alleles.fa"
    fasta.write_text(
        "".join(
            f">{decode_seqid(record.id, default_db='lm').gene}\n{record.sequence}\n"
            for record in iter_fasta(bundled.sequences_path)
        ),
        encoding="utf-8",
    )
    built = runner.invoke(app, ["db", "build", "alleles", str(fasta), "--datadir", str(datadir)])
    assert built.exit_code == 0, built.stderr
    records_path = datadir / "alleles" / "records.jsonl"
    records = [json.loads(line) for line in records_path.read_text(encoding="utf-8").splitlines()]
    merged = [
        {**record, "gene": "merged"} if record["gene"] in {"lmo0737", "lmo1118"} else record
        for record in records
    ]
    records_path.write_text(
        "".join(json.dumps(record) + "\n" for record in merged), encoding="utf-8"
    )

    rebuilt = runner.invoke(app, ["db", "build", "alleles", "--datadir", str(datadir)])

    assert rebuilt.exit_code == 0, rebuilt.stderr
    assert json.loads(rebuilt.stdout)["records"] == 5
    sample = datadir.parent / "sample.fa"
    sample.write_text(
        "".join(
            f">{record.id}\n{record.sequence}\n"
            for record in iter_fasta(datadir / "alleles" / "sequences")
            if decode_seqid(record.id, default_db="alleles").gene == "merged"
        ),
        encoding="utf-8",
    )
    screened = runner.invoke(
        app,
        ["screen", str(sample), "--db", "alleles", "--datadir", str(datadir), "--quiet"],
    )
    assert screened.exit_code == 0, screened.stderr
    genes = {line.split("\t")[5] for line in screened.stdout.splitlines()[1:]}
    assert genes == {"merged"}


def test_rebuild_rejects_a_cluster_database(datadir: Path) -> None:
    """Given a cluster database, When `db build NAME` (no INPUT) runs, Then exit 4
    DB_REBUILD_CLUSTER points at `db build NAME FILE` from the source GBK/GFF."""
    materialize_bundled("vp_ok", datadir, quiet=True)

    result = runner.invoke(app, ["db", "build", "vp_ok", "--datadir", str(datadir)])

    assert result.exit_code == 4
    envelope = json.loads(result.stderr.splitlines()[-1])
    assert envelope["code"] == "DB_REBUILD_CLUSTER"


def test_rebuild_rejects_an_unknown_database(datadir: Path) -> None:
    """Given a name that is not installed, Then exit 4 DATABASE_NOT_FOUND."""
    result = runner.invoke(app, ["db", "build", "nope", "--datadir", str(datadir)])
    assert result.exit_code == 4
    assert json.loads(result.stderr.splitlines()[-1])["code"] == "DATABASE_NOT_FOUND"


def test_a_broken_record_edit_leaves_the_database_intact(datadir: Path) -> None:
    """Given a records file with an unparseable line, When the in-place rebuild runs,
    Then it fails without touching the installed artifacts (staged publish)."""
    _materialize("lm_doumith", "lm", datadir)
    before = (datadir / "lm" / "sequences").read_bytes()
    records_path = datadir / "lm" / "records.jsonl"
    records_path.write_text(
        records_path.read_text(encoding="utf-8") + "{not json\n", encoding="utf-8"
    )

    result = runner.invoke(app, ["db", "build", "lm", "--datadir", str(datadir)])

    assert result.exit_code == 5
    assert json.loads(result.stderr.splitlines()[-1])["code"] == "RECORDS_MALFORMED"
    assert (datadir / "lm" / "sequences").read_bytes() == before


def test_rebuild_preserves_the_info_sidecar(datadir: Path) -> None:
    """Given a database carrying an info.json sidecar, When it is rebuilt,
    Then the sidecar survives — a curation edit must not drop the database's
    own literature and caveats."""
    _materialize("lm_doumith", "lm", datadir)
    sidecar = datadir / "lm" / "info.json"
    sidecar.write_text(
        json.dumps({"references": [{"citation": "Smith 2026"}], "caveats": ["keep me"]}),
        encoding="utf-8",
    )

    result = runner.invoke(app, ["db", "build", "lm", "--datadir", str(datadir)])

    assert result.exit_code == 0, result.stderr
    assert json.loads(sidecar.read_text(encoding="utf-8"))["caveats"] == ["keep me"]


def test_in_place_build_rejects_input_shaping_flags(datadir: Path) -> None:
    """Given no INPUT, When a new-database flag is passed, Then the call fails
    as a usage error naming the flags instead of silently ignoring them."""
    _materialize("lm_doumith", "lm", datadir)
    tsv = datadir.parent / "meta.tsv"
    tsv.write_text("gene\taccession\n", encoding="utf-8")

    result = runner.invoke(app, ["db", "build", "lm", "--tsv", str(tsv), "--datadir", str(datadir)])

    assert result.exit_code == 2
    envelope = json.loads(result.stderr.splitlines()[-1])
    assert envelope["code"] == "USAGE_ERROR"
    assert "--tsv" in envelope["message"]
