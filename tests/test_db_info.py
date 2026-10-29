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
from gapit.db_info_manual import render_db_manual
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
    result = runner.invoke(app, ["db", "info", "-d", "cje_cps", "--datadir", str(datadir)])

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

    as_json = runner.invoke(
        app, ["db", "info", "-d", "cje_cps", "--datadir", str(datadir), "--json"]
    )
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
    assert document.records == 42
    assert document.kind == "gene"
    assert document.typing_schema == "gapit.typing/2"
    assert document.build is not None and document.build.tool.startswith("gapit ")
    assert document.snapshot == "2026-10-29"
    assert len(document.references) == 4
    assert len(document.caveats) >= 5
    assert document.destination == str(datadir / "cje_cps")


def test_info_for_a_catalog_provider_is_available_and_licensed(datadir: Path) -> None:
    """Given a registry provider that is not installed, When `db info` runs,
    Then the row reads `available`, carries the catalog licence, and reports no
    record count or build provenance."""
    result = runner.invoke(app, ["db", "info", "-d", "card", "--datadir", str(datadir), "--json"])

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

    result = runner.invoke(app, ["db", "info", "-d", "labgenes", "--datadir", str(datadir)])

    assert result.exit_code == 4
    envelope = json.loads(result.stderr.splitlines()[-1])
    assert envelope["code"] == "DBINFO_MALFORMED"


def test_info_rejects_an_unknown_database(datadir: Path) -> None:
    """Given a name that is neither catalog, bundled, nor installed, Then exit
    4 DATABASE_NOT_FOUND pointing at `db list`."""
    result = runner.invoke(app, ["db", "info", "-d", "nope", "--datadir", str(datadir)])

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


def test_manual_rendering_is_a_man_page(datadir: Path) -> None:
    """Given an installed panel, When the man-style rendering is requested,
    Then the page carries the classic header, the uppercase sections in order,
    a hanging-indent numbered reference, and no ANSI escapes."""
    materialize_bundled("cje_cps", datadir, quiet=True)
    document = collect_db_info("cje_cps", datadir)

    page = render_db_manual(document, width=80)

    lines = page.splitlines()
    assert lines[0].startswith("CJE_CPS(1)") and lines[0].endswith("CJE_CPS(1)")
    assert "gapit database manual" in lines[0]
    headings = [line for line in lines if line and not line.startswith(" ")]
    assert headings == [
        "CJE_CPS(1)                   gapit database manual                    CJE_CPS(1)",
        "NAME",
        "SYNOPSIS",
        "DESCRIPTION",
        "STATUS",
        "PROVENANCE",
        "REFERENCES",
        "CAVEATS",
        "USAGE",
        "SEE ALSO",
    ]
    body = [line for line in lines if line.startswith(" " * 7)]
    assert any(line.startswith("       cje_cps — Campylobacter") for line in body)
    # Prose sections wrap at the requested width; only the unbreakable
    # filesystem path is allowed to run long.
    prose: list[str] = []
    in_prose = False
    for line in lines:
        if line in {"REFERENCES", "CAVEATS"}:
            in_prose = True
            continue
        if line in {"USAGE", "SEE ALSO"}:
            in_prose = False
        if in_prose and line.startswith(" "):
            prose.append(line)
    assert prose and all(len(line) <= 80 for line in prose)
    assert "       [1] Poly F" in page
    assert "           doi: 10.1371/journal.pone.0144349" in page
    assert "\x1b" not in page


def test_manual_rendering_omits_sections_it_does_not_have(datadir: Path) -> None:
    """Given a bare `db build` product (no manifest prose, no sidecar), When
    rendered, Then only the sections it can fill appear: USAGE is synthesized
    from the kind, while REFERENCES/CAVEATS stay absent instead of printing
    empty headings."""
    fasta = datadir.parent / "genes.fa"
    fasta.write_text(GENE_FASTA, encoding="utf-8")
    built = runner.invoke(app, ["db", "build", "labgenes", str(fasta), "--datadir", str(datadir)])
    assert built.exit_code == 0, built.stderr
    document = collect_db_info("labgenes", datadir)

    page = render_db_manual(document, width=80)

    assert "REFERENCES" not in page
    assert "CAVEATS" not in page
    assert "STATUS" in page and "installed" in page
    # The usage recipe is derived, so even a bare database is runnable.
    assert "USAGE" in page
    assert "$ gapit screen -d labgenes [FILE]... --format json" in page
    # No build-less provenance rows are invented.
    assert "Built by" in page  # this one IS installed, so it has build provenance
    assert "SEE ALSO" in page


def test_usage_lines_are_not_broken_mid_command(datadir: Path) -> None:
    """Given a usage recipe with trailing comments, When rendered, Then each
    command stays on one unwrapped line and its comment moves below."""
    materialize_bundled("cje_cps", datadir, quiet=True)
    document = collect_db_info("cje_cps", datadir)

    page = render_db_manual(document, width=80)

    assert "       $ gapit screen isolate.fa --db cje_cps --output cps.tsv\n" in page
    assert "         # stage 1: detect CPS amplicons" in page


def test_tsv_flag_forces_the_table_on_a_terminal(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a TTY stdout, When --tsv is passed, Then the stable FIELD/VALUE
    table is printed instead of the man page (agents and diffs need it)."""
    monkeypatch.setattr("gapit.cmd_db_info._stdout_is_tty", lambda: True)

    forced = runner.invoke(app, ["db", "info", "-d", "cje_cps", "--datadir", str(datadir), "--tsv"])

    assert forced.exit_code == 0, forced.stderr
    assert forced.stdout.startswith("NAME\tcje_cps\n")
    assert "SYNOPSIS" not in forced.stdout


def test_tty_stdout_gets_the_manual(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given a TTY stdout and no --tsv, When `db info` runs, Then the manual is
    printed (the terminal default)."""
    monkeypatch.setattr("gapit.cmd_db_info._stdout_is_tty", lambda: True)

    result = runner.invoke(app, ["db", "info", "-d", "cje_cps", "--datadir", str(datadir)])

    assert result.exit_code == 0, result.stderr
    assert "CJE_CPS(1)" in result.stdout
    assert "NAME\n       cje_cps —" in result.stdout


def test_every_catalog_and_bundled_database_carries_guidance(datadir: Path) -> None:
    """Given the whole catalog, When each database's info is collected, Then
    every registry provider and every bundled panel reports at least one
    literature reference and one interpretation caveat, and a synthesized
    usage recipe — the completeness gate that stops a database shipping
    without documentation."""
    from gapit.bundled import bundled_names
    from gapit.providers import REGISTRY

    missing: list[str] = []
    for name in sorted(set(REGISTRY) | set(bundled_names())):
        document = collect_db_info(name, datadir)
        if not document.references or not document.caveats or not document.usage:
            missing.append(
                f"{name} (refs={len(document.references)}, "
                f"caveats={len(document.caveats)}, usage={len(document.usage)})"
            )
    assert missing == [], f"databases without db info guidance: {missing}"


def test_every_database_has_at_least_one_resolvable_reference(datadir: Path) -> None:
    """Given each catalog/bundled database, When its references are inspected,
    Then at least one carries a DOI or URL an agent can follow. Individual
    entries may be bare citations (a national standard, an older paper without
    a pinned DOI), but a page never points nowhere."""
    from gapit.bundled import bundled_names
    from gapit.providers import REGISTRY

    unresolved: list[str] = []
    for name in sorted(set(REGISTRY) | set(bundled_names())):
        references = collect_db_info(name, datadir).references
        if references and not any(reference.doi or reference.url for reference in references):
            unresolved.append(name)
    assert unresolved == [], f"databases whose references are all unfollowable: {unresolved}"


def test_bundled_cje_cps_usage_is_declared_not_synthesized(datadir: Path) -> None:
    """Given cje_cps, When info is collected, Then its declared two-stage
    recipe (screen to a table, then typing) is reported verbatim rather than
    replaced by the generic kind default."""
    document = collect_db_info("cje_cps", datadir)

    assert any("--output cps.tsv" in step for step in document.usage)
    assert any(step.startswith("gapit typing") for step in document.usage)
