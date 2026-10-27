"""`gapit db info`: one database's provenance, references, usage, and caveats.

The command is a pure read over four existing sources (provider catalog,
``bundled.json``, the installed manifest, and an optional ``info.json``
sidecar), so these tests pin the merge order rather than the content: a
bundled-but-unmaterialized panel reports ``bundled``, materializing it flips
the same command to ``installed`` with build provenance, a catalog provider
stays ``available``, and a user-built database can carry its own literature.
"""

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from gapit.bundled import materialize_bundled
from gapit.cli import app
from gapit.db_info_ops import collect_db_info, render_db_info

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"
GENE_FASTA = ">labcur1 cur1\nATGGCTAGCTAGGCTAAGGCTTAGCGGATCCGATTACAGGCATCAGGGCTTAACC\n"

runner = CliRunner()


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    target = tmp_path / "datadir"
    target.mkdir()
    return target


def test_info_for_a_bundled_database_carries_literature_and_caveats(
    datadir: Path,
) -> None:
    """Given the bundled cje_cps database (not yet materialized), When
    `db info` runs, Then the status is `bundled`, the paper reference and the
    measured caveats are reported, and the JSON document validates against
    the registered schema."""
    result = runner.invoke(app, ["db", "info", "cje_cps", "--datadir", str(datadir)])

    assert result.exit_code == 0, result.stderr
    rows = [line.split("\t") for line in result.stdout.splitlines()]
    fields = {field for field, _value in rows}
    assert {
        "NAME",
        "DESCRIPTION",
        "KIND",
        "DBTYPE",
        "STATUS",
        "REFERENCE",
        "CAVEAT",
        "USAGE",
    } <= fields
    status = next(value for field, value in rows if field == "STATUS")
    assert status == "bundled"
    references = [value for field, value in rows if field == "REFERENCE"]
    assert any("Poly" in reference and "e0144349" in reference for reference in references)
    assert any("Klena" in reference and "lpxA" in reference for reference in references)
    caveats = [value for field, value in rows if field == "CAVEAT"]
    assert any("47 Penner serotypes collapse into 35 capsule types" in caveat for caveat in caveats)
    assert any("E-value bound" in caveat for caveat in caveats)
    # A reference's DOI/note rows follow its own REFERENCE row.
    first_reference = rows.index(next(row for row in rows if row[0] == "REFERENCE"))
    assert rows[first_reference + 1][0] == "REFERENCE_DOI"

    as_json = runner.invoke(app, ["db", "info", "cje_cps", "--datadir", str(datadir), "--json"])
    assert as_json.exit_code == 0, as_json.stderr
    document = json.loads(as_json.stdout)
    assert document["schema"] == "gapit.dbinfo/1"
    assert document["name"] == "cje_cps"
    assert document["status"] == "bundled"
    assert len(document["references"]) == 4
    assert document["references"][0]["doi"] == "10.1371/journal.pone.0144349"
    assert "records" not in document
    assert "build" not in document


def test_info_flips_to_installed_with_build_provenance(datadir: Path) -> None:
    """Given the same bundled panel materialized into the datadir, When
    `db info` runs, Then the status and record count come from the manifest,
    the provenance block appears, and the bundled guidance survives."""
    materialize_bundled("cje_cps", datadir, quiet=True)

    document = collect_db_info("cje_cps", datadir)

    assert document.status == "installed"
    assert document.records == 37
    assert document.kind == "gene"
    assert document.typing_schema == "gapit.typing/2"
    assert document.build is not None and document.build.tool.startswith("gapit ")
    assert document.snapshot == "2026-10-27"
    assert len(document.references) == 4
    assert len(document.caveats) >= 5
    assert document.destination == str(datadir / "cje_cps")


def test_info_for_a_catalog_provider_is_available_and_licensed(datadir: Path) -> None:
    """Given a registry provider that is not installed, When `db info` runs,
    Then the row reads `available`, carries the catalog licence, and reports no
    record count or build provenance."""
    result = runner.invoke(app, ["db", "info", "card", "--datadir", str(datadir), "--json"])

    assert result.exit_code == 0, result.stderr
    document = json.loads(result.stdout)
    assert document["status"] == "available"
    assert document["vendor"] == "McMaster University"
    assert "non-commercial" in document["license"]
    assert "records" not in document


def test_info_reads_a_user_database_sidecar(datadir: Path) -> None:
    """Given a `db build` product with an info.json sidecar, When `db info`
    runs, Then the sidecar's references and caveats are reported."""
    fasta = datadir.parent / "genes.fa"
    fasta.write_text(GENE_FASTA, encoding="utf-8")
    built = runner.invoke(app, ["db", "build", "labgenes", str(fasta), "--datadir", str(datadir)])
    assert built.exit_code == 0, built.stderr
    sidecar = datadir / "labgenes" / "info.json"
    sidecar.write_text(
        json.dumps(
            {
                "references": [{"citation": "Smith 2026, Lab panel (unpublished)"}],
                "caveats": ["unpublished panel: validate before reporting"],
                "usage": ["gapit screen x.fa --db labgenes"],
            }
        ),
        encoding="utf-8",
    )

    document = collect_db_info("labgenes", datadir)

    assert document.status == "installed"
    assert document.kind == "gene"
    assert document.source_urls == ("local",)
    assert [reference.citation for reference in document.references] == [
        "Smith 2026, Lab panel (unpublished)"
    ]
    assert document.caveats == ("unpublished panel: validate before reporting",)
    assert document.usage == ("gapit screen x.fa --db labgenes",)


def test_info_rejects_a_malformed_sidecar(datadir: Path) -> None:
    """Given an unparseable info.json, When `db info` runs, Then exit 4
    DBINFO_MALFORMED names the file instead of silently dropping provenance."""
    fasta = datadir.parent / "genes.fa"
    fasta.write_text(GENE_FASTA, encoding="utf-8")
    built = runner.invoke(app, ["db", "build", "labgenes", str(fasta), "--datadir", str(datadir)])
    assert built.exit_code == 0, built.stderr
    (datadir / "labgenes" / "info.json").write_text("{not json", encoding="utf-8")

    result = runner.invoke(app, ["db", "info", "labgenes", "--datadir", str(datadir)])

    assert result.exit_code == 4
    envelope = json.loads(result.stderr.splitlines()[-1])
    assert envelope["code"] == "DBINFO_MALFORMED"


def test_info_rejects_an_unknown_database(datadir: Path) -> None:
    """Given a name that is neither catalog, bundled, nor installed, Then exit
    4 DATABASE_NOT_FOUND pointing at `db list`."""
    result = runner.invoke(app, ["db", "info", "nope", "--datadir", str(datadir)])

    assert result.exit_code == 4
    envelope = json.loads(result.stderr.splitlines()[-1])
    assert envelope["code"] == "DATABASE_NOT_FOUND"


def test_info_render_is_two_columns_and_tab_safe(datadir: Path) -> None:
    """Given a document with multi-line caveat text, When rendered, Then every
    row is a two-cell TSV line (tabs and newlines in values are sanitized)."""
    materialize_bundled("cje_cps", datadir, quiet=True)
    document = collect_db_info("cje_cps", datadir).model_copy(
        update={"caveats": ("two\tcells\nand a newline",)}
    )

    rendered = render_db_info(document)

    rows = rendered.splitlines()
    assert all(len(row.split("\t")) == 2 for row in rows)
    assert "\nand a newline" not in rendered


def test_dbinfo_schema_is_registered() -> None:
    """Given `gapit schema dbinfo`, When run, Then the JSON schema of
    gapit.dbinfo/1 is printed — the document is part of the public surface."""
    result = runner.invoke(app, ["schema", "dbinfo"])

    assert result.exit_code == 0, result.stderr
    schema = json.loads(result.stdout)
    assert schema["title"] == "DbInfoDocument"
