"""`gspec typing` on cluster screen JSON (already-evaluated phenotypes)."""

from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from gspec.cli import app

runner = CliRunner()


def _call(phenotype: str, score: float = 1.0) -> dict[str, object]:
    return {
        "phenotype": phenotype,
        "score": score,
        "confidence": "high",
        "components": [],
    }


def _document(
    schema: str = "gapit.cluster/2", *, db: str = "vpa_ok", file: str = "sample.fna"
) -> dict[str, object]:
    return {
        "schema": schema,
        "params": {"db": db},
        "files": [
            {
                "file": file,
                "phenotypes": {
                    "o_group": _call("O3"),
                    "k_group": _call("K6"),
                    "serotype": _call("O3:K6"),
                },
            }
        ],
    }


def _write(tmp_path: Path, document: dict[str, object], name: str = "cluster.json") -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(document), encoding="utf-8")
    return path


def test_typing_reads_cluster2_json(tmp_path: Path) -> None:
    """Given a gapit.cluster/2 JSON, When typed, Then the concise table is
    rendered from its already-evaluated phenotypes."""
    path = _write(tmp_path, _document())
    result = runner.invoke(app, ["typing", str(path), "--quiet"])
    assert result.exit_code == 0, result.stderr
    lines = result.stdout.splitlines()
    assert lines[0].split("\t")[:3] == ["FILE", "SCHEME", "PHENOTYPE"]
    assert [line.split("\t")[2] for line in lines[1:]] == ["O3", "K6", "O3:K6"]


def test_typing_reads_locus_map_json(tmp_path: Path) -> None:
    """Given a gapit.locus_map/1 document, When typed, Then its phenotypes
    map is rendered without needing the database or a second screen."""
    path = _write(tmp_path, _document("gapit.locus_map/1"), "detail.json")
    result = runner.invoke(app, ["typing", str(path), "--format", "json", "--quiet"])
    assert result.exit_code == 0, result.stderr
    document = json.loads(result.stdout)
    assert document["schema"] == "gapit.typing_result/1"
    assert document["db"] == "vpa_ok"
    assert document["files"][0]["file"] == "sample.fna"
    assert document["files"][0]["phenotypes"]["serotype"]["phenotype"] == "O3:K6"


def test_typing_rejects_mixed_cluster_databases(tmp_path: Path) -> None:
    """Given two cluster documents from different databases, When typed,
    Then fail with DATABASE_MISMATCH instead of silently merging."""
    first = _write(tmp_path, _document(db="vpa_ok"), "a.json")
    second = _write(tmp_path, _document(db="other"), "b.json")
    result = runner.invoke(app, ["typing", str(first), str(second), "--quiet"])
    assert result.exit_code == 5
    envelope = json.loads(result.stderr.splitlines()[-1])
    assert envelope["code"] == "DATABASE_MISMATCH"


def test_typing_rejects_unknown_json_schema(tmp_path: Path) -> None:
    """Given a JSON document that is not a cluster result, When typed, Then
    fail with a typed malformed-input error."""
    path = _write(tmp_path, {"schema": "other/1", "params": {"db": "x"}, "files": []})
    result = runner.invoke(app, ["typing", str(path), "--quiet"])
    assert result.exit_code == 5
    envelope = json.loads(result.stderr.splitlines()[-1])
    assert envelope["code"] == "TYPING_MALFORMED"


def test_typing_rejects_cluster_tsv_with_pointer_to_json(tmp_path: Path) -> None:
    """Given a cluster/2 TSV, When typed, Then fail with a targeted error
    that points at the JSON result carrying scheme metadata."""
    path = tmp_path / "cluster.tsv"
    path.write_text("FILE\tSEROTYPE\tO_GROUP\tK_GROUP\n", encoding="utf-8")
    result = runner.invoke(app, ["typing", str(path), "--quiet"])
    assert result.exit_code == 5
    envelope = json.loads(result.stderr.splitlines()[-1])
    assert envelope["code"] == "TYPING_CLUSTER_TABLE"
    assert "format json" in envelope["message"]
