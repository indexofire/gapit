"""Pre-0.9 (``gapit.*``) spellings still load and are re-emitted as ``gspec.*``.

0.9 renamed every on-disk contract: schema ids, the manifest filename, the
FASTA header prefix and the ``header_format`` value. Existing databases,
archived JSON and user-authored sidecars keep working because every reader
accepts the legacy spelling and every writer emits the current one. These
tests pin that promise.
"""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

from pydantic import TypeAdapter

from gspec.dbcodec import PREFIX, decode_seqid, encode_seqid, is_gspec_header
from gspec.formats.json import ReportDocument
from gspec.layout import LEGACY_MANIFEST_NAME, MANIFEST_NAME, manifest_path
from gspec.records import Manifest, read_manifest, write_manifest
from gspec.typing_models import read_typing_document, typing_schema_of

report_adapter = TypeAdapter(ReportDocument)
NOW = datetime(2026, 11, 7, 12, 0, 0, tzinfo=UTC)


def _minimal_manifest(**overrides: object) -> dict[str, object]:
    document: dict[str, object] = {
        "schema": "gapit.manifest/1",
        "name": "panel",
        "source_urls": ["local"],
        "fetched_at": "2026-11-07T12:00:00Z",
        "sha256": "0" * 64,
        "n_records": 1,
        "dbtype": "nucl",
    }
    document.update(overrides)
    return document


def test_legacy_schema_id_is_accepted_and_canonicalised(tmp_path: Path) -> None:
    path = tmp_path / "report.json"
    path.write_text(
        json.dumps(
            {
                "schema": "gapit.report/1",
                "tool": {"name": "gspec", "version": "0.8.1"},
                "created_at": "2026-11-07T12:00:00Z",
                "params": {"db": "panel", "minid": 80.0, "mincov": 80.0, "threads": 1},
                "files": [],
            }
        ),
        encoding="utf-8",
    )
    document = report_adapter.validate_json(path.read_text(encoding="utf-8"))
    assert document.schema_name == "gspec.report/1"
    assert document.model_dump(by_alias=True)["schema"] == "gspec.report/1"


def test_legacy_manifest_filename_is_read(tmp_path: Path) -> None:
    db_dir = tmp_path / "panel"
    db_dir.mkdir()
    legacy = db_dir / LEGACY_MANIFEST_NAME
    legacy.write_text(json.dumps(_minimal_manifest()), encoding="utf-8")
    assert manifest_path(db_dir) == legacy
    assert read_manifest(legacy).name == "panel"


def test_current_manifest_filename_wins(tmp_path: Path) -> None:
    db_dir = tmp_path / "panel"
    db_dir.mkdir()
    (db_dir / LEGACY_MANIFEST_NAME).write_text(json.dumps(_minimal_manifest()), encoding="utf-8")
    (db_dir / MANIFEST_NAME).write_text(
        json.dumps(_minimal_manifest(name="preferred")), encoding="utf-8"
    )
    assert manifest_path(db_dir) == db_dir / MANIFEST_NAME
    assert read_manifest(manifest_path(db_dir)).name == "preferred"


def test_manifest_writes_the_current_filename(tmp_path: Path) -> None:
    manifest = Manifest.model_validate(_minimal_manifest())
    target = tmp_path / MANIFEST_NAME
    write_manifest(manifest, target)
    assert target.read_text(encoding="utf-8").startswith('{\n  "schema": "gspec.manifest/1"')


def test_legacy_header_format_is_accepted(tmp_path: Path) -> None:
    manifest = Manifest.model_validate(_minimal_manifest(header_format="gapit/v1"))
    assert manifest.header_format == "gspec/v1"


def test_header_prefix_round_trips_and_decodes_the_legacy_spelling() -> None:
    assert PREFIX == "gspec|"
    seqid = encode_seqid("panel", "gene1", "ACC", ["function"])
    assert seqid.startswith("gspec|db=panel|gene=gene1|")
    header = decode_seqid(seqid, "fallback")
    assert (header.gene, header.accession) == ("gene1", "ACC")

    legacy = seqid.replace("gspec|", "gapit|", 1)
    assert is_gspec_header(legacy)
    assert decode_seqid(legacy, "fallback").gene == "gene1"


def test_legacy_typing_document_loads_and_is_reported_canonically(tmp_path: Path) -> None:
    path = tmp_path / "typing.json"
    path.write_text(
        json.dumps(
            {
                "schema": "gapit.typing/2",
                "schemes": [
                    {
                        "name": "serogroup",
                        "cutoff": 0.9,
                        "ambiguity_margin": 0.0,
                        "fallback": "NG",
                        "rules": [
                            {
                                "model": "exact_set",
                                "phenotype": "A",
                                "requires": ["wzx"],
                            }
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    document = read_typing_document(path)
    assert document.schema_name == "gspec.typing/2"
    assert typing_schema_of(path) == "gspec.typing/2"
