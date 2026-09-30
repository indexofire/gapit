"""CLI integration tests: `gapit db build` cluster kind (GBK/GFF inputs).

Real makeblastdb from the pixi env, synthetic fixtures from
tests/data/cluster, everything local to tmp_path — the test_cli_db_build.py
conventions. The gene-kind FASTA path must stay byte-compatible (its own
suite, untouched); these tests pin the NEW cluster artifacts: locus FASTA
``sequences``, ``features.json`` (gapit.features/1), manifest
``kind: cluster``, BLAST index, and the optional validated ``typing.json``.
"""

import json
from pathlib import Path

from pydantic import TypeAdapter
from typer.testing import CliRunner, Result

from gapit.cli import app
from gapit.errors import ErrorEnvelope
from gapit.gbfeatures import FeaturesDocument
from gapit.records import read_manifest
from gapit.typing_models import TypingDocument

DATA = Path(__file__).parent / "data" / "cluster"
DB = "myclusters"

runner = CliRunner()
envelope_adapter = TypeAdapter(ErrorEnvelope)
features_adapter = TypeAdapter(FeaturesDocument)


def last_envelope(stderr: str) -> ErrorEnvelope:
    lines = [line for line in stderr.splitlines() if line.strip()]
    assert lines, "expected an error envelope on stderr"
    return envelope_adapter.validate_json(lines[-1])


def build(*extra: str) -> Result:
    return runner.invoke(app, ["db", "build", *extra])


def test_db_build_gbk_creates_cluster_artifacts(tmp_path: Path) -> None:
    """Given a Bakta-style GBK, When `db build NAME file.gbk`, Then exit 0
    with a gene-shaped receipt (records = loci), the locus FASTA carries one
    record per locus, features.json is a valid gapit.features/1, the manifest
    declares kind cluster with plain headers, the BLAST index exists, and no
    records.jsonl or typing.json is written."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build(DB, str(DATA / "bakta_style.gbk"), "--datadir", str(datadir))

    assert result.exit_code == 0, result.stderr
    db_dir = datadir / DB
    assert json.loads(result.stdout) == {
        "db": DB,
        "records": 2,
        "dbtype": "nucl",
        "destination": str(db_dir),
    }
    sequences = (db_dir / "sequences").read_text(encoding="utf-8")
    assert [line for line in sequences.splitlines() if line.startswith(">")] == [
        ">KL9001",
        ">KL9002",
    ]
    document = features_adapter.validate_json(
        (db_dir / "features.json").read_text(encoding="utf-8")
    )
    assert document.schema_name == "gapit.features/1"
    assert [locus.id for locus in document.loci] == ["KL9001", "KL9002"]
    assert [len(locus.genes) for locus in document.loci] == [3, 3]
    manifest = read_manifest(db_dir / "gapit-manifest.json")
    assert manifest.kind == "cluster"
    assert manifest.header_format == "plain"
    assert manifest.n_records == 2
    assert manifest.dbtype == "nucl"
    assert (db_dir / "sequences.nin").is_file()
    assert not (db_dir / "records.jsonl").exists()
    assert not (db_dir / "typing.json").exists()


def test_db_build_gff3_sidecar_creates_cluster_artifacts(tmp_path: Path) -> None:
    """Given a GFF3 whose sequences come from a sidecar .fa, When built,
    Then the cluster artifacts carry the FASTA-defined loci and their CDS
    genes (coordinate semantics identical to the GBK path)."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build(DB, str(DATA / "clusters_sidecar.gff3"), "--datadir", str(datadir))

    assert result.exit_code == 0, result.stderr
    db_dir = datadir / DB
    document = features_adapter.validate_json(
        (db_dir / "features.json").read_text(encoding="utf-8")
    )
    assert [locus.id for locus in document.loci] == ["cl1", "cl2"]
    assert document.loci[1].genes[0].gene_id == "bigA"
    assert read_manifest(db_dir / "gapit-manifest.json").kind == "cluster"
    assert (db_dir / "sequences.nin").is_file()


def test_db_build_gzipped_gbk(tmp_path: Path) -> None:
    """Given a .gz GBK, When built, Then the cluster build decompresses
    transparently."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build(DB, str(DATA / "gz_input.gbk.gz"), "--datadir", str(datadir))

    assert result.exit_code == 0, result.stderr
    assert read_manifest(datadir / DB / "gapit-manifest.json").n_records == 1


def test_db_build_kind_flag_must_match_detected_format(tmp_path: Path) -> None:
    """Given a GBK input and --kind gene (or a FASTA input and --kind
    cluster), When built, Then exit 2 with a usage envelope naming both the
    requested and the detected kind; a matching --kind cluster succeeds."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fasta = tmp_path / "genes.fa"
    fasta.write_text(">g1 demo\nACGTACGTACGTACGTACGTACGTACGTACGT\n", encoding="utf-8")

    contradicted = build(
        DB, str(DATA / "bakta_style.gbk"), "--datadir", str(datadir), "--kind", "gene"
    )
    assert contradicted.exit_code == 2
    envelope = last_envelope(contradicted.stderr)
    assert envelope.code == "USAGE_ERROR"
    assert envelope.context["kind"] == "gene"
    assert envelope.context["detected"] == "cluster"

    reversed_case = build("faclusters", str(fasta), "--datadir", str(datadir), "--kind", "cluster")
    assert reversed_case.exit_code == 2
    assert last_envelope(reversed_case.stderr).code == "USAGE_ERROR"

    matching = build(
        "explicit", str(DATA / "bakta_style.gbk"), "--datadir", str(datadir), "--kind", "cluster"
    )
    assert matching.exit_code == 0, matching.stderr
    assert read_manifest(datadir / "explicit" / "gapit-manifest.json").kind == "cluster"


def test_db_build_gene_kind_stays_gene(tmp_path: Path) -> None:
    """Given a plain FASTA (unknown-to-cluster suffixes), When built, Then
    the manifest declares kind gene with gapit/v1 headers — the gene pipeline
    is the untouched default."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fasta = tmp_path / "genes.fa"
    fasta.write_text(">demov2 demo protein\nACGTACGTACGTACGTACGTACGTACGTACGT\n", encoding="utf-8")

    result = build(DB, str(fasta), "--datadir", str(datadir))

    assert result.exit_code == 0, result.stderr
    manifest = read_manifest(datadir / DB / "gapit-manifest.json")
    assert manifest.kind == "gene"
    assert manifest.header_format == "gapit/v1"
    assert not (datadir / DB / "features.json").exists()


def test_db_build_typing_validated_and_copied(tmp_path: Path) -> None:
    """Given a valid --typing document, When built, Then typing.json lands
    in the db directory byte-identical to the source."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build(
        DB,
        str(DATA / "bakta_style.gbk"),
        "--datadir",
        str(datadir),
        "--typing",
        str(DATA / "typing_valid.json"),
    )

    assert result.exit_code == 0, result.stderr
    copied = datadir / DB / "typing.json"
    assert copied.read_bytes() == (DATA / "typing_valid.json").read_bytes()
    document = TypingDocument.model_validate_json(copied.read_text(encoding="utf-8"))
    assert document.fallback == "unknown"


def test_db_build_typing_invalid_stops_the_build(tmp_path: Path) -> None:
    """Given an invalid --typing document, When built, Then exit 4 with
    TYPING_MALFORMED and no manifest certifies a half-built database."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build(
        DB,
        str(DATA / "bakta_style.gbk"),
        "--datadir",
        str(datadir),
        "--typing",
        str(DATA / "typing_invalid.json"),
    )

    assert result.exit_code == 4
    envelope = last_envelope(result.stderr)
    assert envelope.code == "TYPING_MALFORMED"
    assert not (datadir / DB / "gapit-manifest.json").exists()


def test_db_build_typing_with_gene_kind_is_usage_error(tmp_path: Path) -> None:
    """Given a FASTA gene build with --typing, When built, Then exit 2 —
    the option belongs to cluster builds only."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fasta = tmp_path / "genes.fa"
    fasta.write_text(">g1 demo\nACGTACGTACGTACGTACGTACGTACGTACGT\n", encoding="utf-8")

    result = build(
        DB, str(fasta), "--datadir", str(datadir), "--typing", str(DATA / "typing_valid.json")
    )

    assert result.exit_code == 2
    assert last_envelope(result.stderr).code == "USAGE_ERROR"


def test_db_build_cluster_input_rejects_gene_options(tmp_path: Path) -> None:
    """Given a GBK input with the FASTA-only --tsv/--dbtype/--description
    options, When built, Then exit 2 usage error rather than silently
    ignoring them."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    tsv = tmp_path / "meta.tsv"
    tsv.write_text("gene\tfunction\ng1\tvirulence\n", encoding="utf-8")

    with_tsv = build(
        DB, str(DATA / "bakta_style.gbk"), "--datadir", str(datadir), "--tsv", str(tsv)
    )
    assert with_tsv.exit_code == 2
    assert last_envelope(with_tsv.stderr).code == "USAGE_ERROR"

    with_dbtype = build(
        DB, str(DATA / "bakta_style.gbk"), "--datadir", str(datadir), "--dbtype", "nucl"
    )
    assert with_dbtype.exit_code == 2

    with_description = build(
        DB, str(DATA / "bakta_style.gbk"), "--datadir", str(datadir), "--description", "x"
    )
    assert with_description.exit_code == 2
