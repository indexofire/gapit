"""gapit.screen/1 per-database defaults: parsing, discovery, validation."""

import json
from pathlib import Path

import pytest

from gapit.db import Database
from gapit.errors import DatabaseError, InputError
from gapit.screen_defaults import (
    SCREEN_FILENAME,
    ScreenDefaults,
    read_screen_defaults,
    screen_defaults,
    screen_defaults_for,
)

VALID = json.dumps(
    {
        "schema": "gapit.screen/1",
        "reads": {
            "min_breadth": 90.0,
            "min_identity": 80.0,
            "min_mapq": 20,
            "min_depth": 5.0,
            "min_relative_depth": 0.1,
            "relative_to": "lpxA",
        },
        "contigs": {"minid": 80.0, "mincov": 80.0, "merge_fragments": True},
    }
)


def test_valid_document_parses() -> None:
    document = ScreenDefaults.model_validate_json(VALID)
    assert document.reads.min_mapq == 20
    assert document.reads.min_relative_depth == 0.1
    assert document.contigs.merge_fragments is True


def test_missing_file_is_input_not_found(tmp_path: Path) -> None:
    with pytest.raises(InputError) as exc_info:
        read_screen_defaults(tmp_path / "nope.json")
    assert exc_info.value.code == "INPUT_NOT_FOUND"


@pytest.mark.parametrize(
    "text",
    [
        "{nope",
        '{"schema": "gapit.screen/2", "reads": {}}',
        '{"schema": "gapit.screen/1", "reads": {"min_mapq": -1}}',
        '{"schema": "gapit.screen/1", "reads": {"min_relative_depth": 1.5}}',
        '{"schema": "gapit.screen/1", "contigs": {"minid": 101}}',
    ],
)
def test_malformed_documents_are_typed(tmp_path: Path, text: str) -> None:
    path = tmp_path / SCREEN_FILENAME
    path.write_text(text, encoding="utf-8")
    with pytest.raises(DatabaseError) as exc_info:
        read_screen_defaults(path)
    assert exc_info.value.code == "SCREEN_DEFAULTS_MALFORMED"


def test_loader_reads_installed_sidecar(tmp_path: Path) -> None:
    (tmp_path / "sequences").write_text(">x\nACGT\n", encoding="utf-8")
    (tmp_path / SCREEN_FILENAME).write_text(VALID, encoding="utf-8")
    database = Database(name="db", path=tmp_path, sequences_path=tmp_path / "sequences")
    document = screen_defaults(database)
    assert document is not None
    assert document.reads.min_identity == 80.0


def test_loader_is_none_without_sidecar(tmp_path: Path) -> None:
    (tmp_path / "sequences").write_text(">x\nACGT\n", encoding="utf-8")
    database = Database(name="db", path=tmp_path, sequences_path=tmp_path / "sequences")
    assert screen_defaults(database) is None


def test_bundled_cje_cps_ships_its_defaults(tmp_path: Path) -> None:
    """The bundled panel carries the curated reads/contigs defaults, so a
    first-use screen does not need flags."""
    document = screen_defaults_for("cje_cps", tmp_path)
    assert document is not None
    assert document.reads.min_identity == 80.0
    assert document.reads.min_mapq == 20
    assert document.reads.min_relative_depth == 0.10
    assert document.contigs.minid == 80.0
    assert document.contigs.merge_fragments is True


def test_db_build_copies_screen_sidecar(tmp_path: Path) -> None:
    """Given a gene build with --screen, When it completes, Then the sidecar
    travels into the database directory and loads as gapit.screen/1."""
    from gapit.db_build_ops import perform_build

    fasta = tmp_path / "genes.fa"
    fasta.write_text(">geneA\nACGTACGTACGT\n>lpxA\nACGTACGTACGT\n", encoding="utf-8")
    screen = tmp_path / "screen.json"
    screen.write_text(VALID, encoding="utf-8")
    datadir = tmp_path / "datadir"
    perform_build(
        "tiny",
        fasta,
        None,
        "nucl",
        "",
        datadir,
        False,
        warn=lambda message: None,
        screen=screen,
    )
    installed = datadir / "tiny" / SCREEN_FILENAME
    assert installed.is_file()
    database = Database(
        name="tiny",
        path=datadir / "tiny",
        sequences_path=datadir / "tiny" / "sequences",
    )
    document = screen_defaults(database)
    assert document is not None
    assert document.contigs.merge_fragments is True


def test_unknown_relative_control_is_a_build_error(tmp_path: Path) -> None:
    from gapit.screen_defaults import validate_screen_defaults

    document = ScreenDefaults.model_validate_json(VALID)
    with pytest.raises(DatabaseError) as exc_info:
        validate_screen_defaults(document, frozenset({"geneA"}), source=tmp_path / "screen.json")
    assert exc_info.value.code == "SCREEN_DEFAULTS_UNKNOWN_GENE"
