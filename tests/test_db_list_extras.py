"""`gapit db list` visibility for databases the registry does not know.

Databases installed into the datadir outside the provider catalog — `db
build` products (gene + cluster kinds) and manifest-less abricate-style
directories — must appear as ``local`` rows in the NAME-sorted listing.
With no extras the output carries no trace of the extension (additive
proof). Real binaries through the CLI surface (test_cli_db_build.py
conventions).
"""

import json
import random
from pathlib import Path

from click.testing import CliRunner, Result
from rich.console import Console

from gapit.cli import app
from gapit.db_ops import db_list_entries, db_list_table
from gapit.providers import REGISTRY
from gapit.records import Manifest, write_manifest

DATA = Path(__file__).parent / "data" / "cluster"

REGISTRY_NAMES = frozenset(REGISTRY)

# Distinct seeds share no 11-mer in practice (test_cli_db_build.py note).
SEQ_A = "".join(random.Random(11).choice("ACGT") for _ in range(240))
SEQ_B = "".join(random.Random(12).choice("ACGT") for _ in range(240))
SEQ_PROT = "".join(random.Random(13).choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(80))

GENE_FASTA = f">demov2 demo beta-lactamase variant 2\n{SEQ_A}\n>demov3 demo efflux pump\n{SEQ_B}\n"

runner = CliRunner()


def build(*extra: str) -> Result:
    return runner.invoke(app, ["db", "build", *extra])


def list_rows(datadir: Path) -> list[str]:
    result = runner.invoke(app, ["db", "list", "--datadir", str(datadir)])
    assert result.exit_code == 0, result.stderr
    return result.stdout.splitlines()


def test_db_list_includes_built_gene_and_cluster_databases(tmp_path: Path) -> None:
    """Given a datadir holding one gene-kind and one cluster-kind `db build`
    product, When `db list`, Then both render as local rows in the globally
    NAME-sorted listing — STATUS from the manifest record count, PROVIDER
    local, DBTYPE nucl, empty DESCRIPTION (no manifest note)."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fasta = tmp_path / "genes.fa"
    fasta.write_text(GENE_FASTA, encoding="utf-8")
    assert build("myamr", str(fasta), "--datadir", str(datadir)).exit_code == 0
    assert build("myloci", str(DATA / "bakta_style.gbk"), "--datadir", str(datadir)).exit_code == 0

    lines = list_rows(datadir)

    assert lines[0] == "NAME\tDESCRIPTION\tSTATUS\tDBTYPE\tTYPING\tPROVIDER"
    rows = [line.split("\t") for line in lines[1:]]
    names = [row[0] for row in rows]
    assert names == sorted(names)
    by_name = {row[0]: row for row in rows}
    assert by_name["myamr"] == ["myamr", "", "installed (2)", "nucl", "-", "local"]
    assert by_name["myloci"] == ["myloci", "", "installed (2)", "nucl", "-", "local"]


def test_db_list_json_marks_extras_local_and_keeps_registry_entries_untouched(
    tmp_path: Path,
) -> None:
    """Given the same mixed datadir, When `db list --json`, Then the extra
    entries sit inside providers with the additive source=local field (plus
    vendor local and the manifest kind), while every registry entry simply
    lacks the field — the absent-means-registry rule (the wheel's bundled
    ecoli_dec row carries source=bundled and is excluded from that check)."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fasta = tmp_path / "genes.fa"
    fasta.write_text(GENE_FASTA, encoding="utf-8")
    assert build("myamr", str(fasta), "--datadir", str(datadir)).exit_code == 0
    assert build("myloci", str(DATA / "bakta_style.gbk"), "--datadir", str(datadir)).exit_code == 0

    result = runner.invoke(app, ["db", "list", "--datadir", str(datadir), "--json"])
    assert result.exit_code == 0
    document = json.loads(result.stdout)
    assert document["schema"] == "gapit.dblist/1"
    providers = document["providers"]
    by_name = {entry["name"]: entry for entry in providers}
    assert by_name["myamr"]["source"] == "local"
    assert by_name["myamr"]["vendor"] == "local"
    assert by_name["myamr"]["installed"] is True
    assert by_name["myamr"]["records"] == 2
    assert by_name["myamr"]["kind"] == "gene"
    assert by_name["myamr"]["typing"] is False
    assert by_name["myloci"]["source"] == "local"
    assert by_name["myloci"]["kind"] == "cluster"
    assert by_name["myloci"]["typing"] is False
    registry_entries = [entry for entry in providers if "source" not in entry]
    assert registry_entries, "expected the real registry providers"
    assert all(entry["name"] in sorted(REGISTRY_NAMES) for entry in registry_entries)


def test_db_list_json_without_extras_carries_only_the_bundled_source_field(
    tmp_path: Path,
) -> None:
    """Given an empty datadir (no local extras), When `db list --json`, Then
    the only source fields in the document are the eight bundled rows
    (alphabetical, source=bundled, not installed, vendor from bundled.json)
    — registry entries stay field-free and the extension remains additive."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = runner.invoke(app, ["db", "list", "--datadir", str(datadir), "--json"])

    assert result.exit_code == 0
    providers = json.loads(result.stdout)["providers"]
    with_source = [entry for entry in providers if "source" in entry]
    assert [(entry["name"], entry["source"]) for entry in with_source] == [
        ("cje_cps", "bundled"),
        ("ecoh", "bundled"),
        ("ecoli_dec", "bundled"),
        ("lmo_sero", "bundled"),
        ("ncbi", "bundled"),
        ("resfinder", "bundled"),
        ("vpa_ok", "bundled"),
    ]
    for entry in with_source:
        assert entry["installed"] is False
        assert "records" not in entry
        assert entry["dbtype"] == "nucl"
    assert (
        next(entry for entry in with_source if entry["name"] == "ecoli_dec")["vendor"]
        == "gapit-curated (public-domain sources)"
    )


def test_db_list_table_renders_local_rows(tmp_path: Path) -> None:
    """Given a datadir with one built database, When db_list_table renders
    the real entries at width 100, Then the local row appears with its
    STATUS text and the capture carries no ANSI escapes."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fasta = tmp_path / "genes.fa"
    fasta.write_text(GENE_FASTA, encoding="utf-8")
    assert build("myamr", str(fasta), "--datadir", str(datadir)).exit_code == 0

    console = Console(record=True, width=100)
    console.print(db_list_table(db_list_entries(datadir)))
    text = console.export_text()

    assert "myamr" in text
    assert "local" in text
    assert "installed (2)" in text
    assert "\x1b" not in text


def test_db_list_manifestless_dir_counts_fasta_and_detects_dbtype(tmp_path: Path) -> None:
    """Given a datadir with manifest-less directories (abricate-style /
    checksum-installed sequences files), When `db list`, Then each renders
    with the FASTA record count and the dbtype detected from the BLAST index
    suffix — or the mol_type heuristic when no index exists (protein-looking
    letters → prot)."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    legacy = datadir / "alegacy"
    legacy.mkdir()
    (legacy / "sequences").write_text(GENE_FASTA, encoding="utf-8")
    (legacy / "sequences.nin").touch()
    proto = datadir / "zproto"
    proto.mkdir()
    (proto / "sequences").write_text(f">prot_gene\n{SEQ_PROT}\n", encoding="utf-8")

    lines = list_rows(datadir)
    rows = {line.split("\t")[0]: line for line in lines[1:]}

    assert rows["alegacy"] == "alegacy\t\tinstalled (2)\tnucl\t-\tlocal"
    assert rows["zproto"] == "zproto\t\tinstalled (1)\tprot\t-\tlocal"
    names = [line.split("\t")[0] for line in lines[1:]]
    assert names == sorted(names)


def test_db_list_json_surfaces_a_bundled_content_license(tmp_path: Path) -> None:
    """Given the bundled vpa_ok row, When listed, Then its pinned MIT content
    license is carried in the gapit.dblist/1 entry."""
    from gapit.db_ops import db_list_entries

    entries = {entry.name: entry for entry in db_list_entries(tmp_path)}
    assert entries["vpa_ok"].license == "MIT (database content)"
    assert entries["vpa_ok"].typing is True
    assert entries["vpa_ok"].kind == "cluster"
    assert entries["ecoh"].typing is False


def test_db_list_flags_a_locally_built_typed_gene_database(tmp_path: Path) -> None:
    """Given a `db build --typing` gene database, When listed, Then the
    typing capability is true and the TYPING cell names the two-stage
    command path; an untyped built database stays false/`-`."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    typing_data = Path(__file__).parent / "data" / "typing"
    built = build(
        "typed",
        str(typing_data / "markers.fa"),
        "--datadir",
        str(datadir),
        "--typing",
        str(typing_data / "typing_markers_v2.json"),
    )
    assert built.exit_code == 0, built.stderr

    result = runner.invoke(app, ["db", "list", "--datadir", str(datadir), "--json"])
    assert result.exit_code == 0, result.stderr
    entry = next(item for item in json.loads(result.stdout)["providers"] if item["name"] == "typed")
    assert entry["kind"] == "gene"
    assert entry["typing"] is True
    assert "typed\t\tinstalled (4)\tnucl\tcommand\tlocal" in list_rows(datadir)


def test_db_list_does_not_call_a_manifest_only_dir_installed(tmp_path: Path) -> None:
    """Given a manifest without sequences, When db list runs, Then the row is
    available, not installed (a half-built DB is not usable)."""
    datadir = tmp_path / "datadir"
    db_dir = datadir / "card"
    db_dir.mkdir(parents=True)
    write_manifest(
        Manifest(
            name="card",
            source_urls=(),
            fetched_at="2026-01-01T00:00:00Z",
            sha256="0" * 64,
            n_records=1,
            dbtype="nucl",
        ),
        db_dir / "gapit-manifest.json",
    )

    result = runner.invoke(app, ["db", "list", "--datadir", str(datadir), "--json"])

    assert result.exit_code == 0
    by_name = {entry["name"]: entry for entry in json.loads(result.stdout)["providers"]}
    assert by_name["card"]["installed"] is False
    assert "records" not in by_name["card"]


def test_db_list_reports_installed_but_unindexed(tmp_path: Path) -> None:
    """Given a manifest with sequences but no BLAST index, When db list runs,
    Then installed is true but indexed is false (screen still needs db setup)."""
    datadir = tmp_path / "datadir"
    db_dir = datadir / "card"
    db_dir.mkdir(parents=True)
    (db_dir / "sequences").write_text(">g\nACGT\n", encoding="utf-8")
    write_manifest(
        Manifest(
            name="card",
            source_urls=(),
            fetched_at="2026-01-01T00:00:00Z",
            sha256="0" * 64,
            n_records=1,
            dbtype="nucl",
        ),
        db_dir / "gapit-manifest.json",
    )

    result = runner.invoke(app, ["db", "list", "--datadir", str(datadir), "--json"])

    assert result.exit_code == 0
    by_name = {entry["name"]: entry for entry in json.loads(result.stdout)["providers"]}
    assert by_name["card"]["installed"] is True
    assert by_name["card"]["indexed"] is False


def test_db_list_name_prints_space_separated_sorted_names(tmp_path: Path) -> None:
    """Given a catalog, When db list -n runs, Then stdout is one space-separated
    line of names in NAME order, with no header and no TSV fields."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    entries = db_list_entries(datadir)
    names = [entry.name for entry in entries]

    result = runner.invoke(app, ["db", "list", "-n", "--datadir", str(datadir)])

    assert result.exit_code == 0
    assert result.stdout == " ".join(names) + "\n"
    assert names == sorted(names)
    assert "\t" not in result.stdout
    assert "NAME" not in result.stdout


def test_db_list_name_and_json_are_mutually_exclusive(tmp_path: Path) -> None:
    """Given -n and --json together, When db list runs, Then it is a usage
    error instead of silently choosing one rendering."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    result = runner.invoke(app, ["db", "list", "-n", "--json", "--datadir", str(datadir)])
    assert result.exit_code == 2
    assert "not both" in result.stderr


def test_db_list_reads_local_info_description_and_vendor(tmp_path: Path) -> None:
    """Given a local database with an info.json description/vendor, When
    listed, Then the row carries both instead of the anonymous local label."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fasta = tmp_path / "genes.fa"
    fasta.write_text(GENE_FASTA, encoding="utf-8")
    assert build("myamr", str(fasta), "--datadir", str(datadir)).exit_code == 0
    (datadir / "myamr" / "info.json").write_text(
        json.dumps(
            {
                "description": "Local AMR panel",
                "vendor": "gapit(Hangzhou CDC)",
            }
        ),
        encoding="utf-8",
    )

    result = runner.invoke(app, ["db", "list", "--datadir", str(datadir), "--json"])
    assert result.exit_code == 0, result.stderr
    entry = next(item for item in json.loads(result.stdout)["providers"] if item["name"] == "myamr")
    assert entry["description"] == "Local AMR panel"
    assert entry["vendor"] == "gapit(Hangzhou CDC)"

    rows = list_rows(datadir)
    assert "myamr\tLocal AMR panel\tinstalled (2)\tnucl\t-\tgapit(Hangzhou CDC)" in rows


def test_db_list_malformed_local_info_is_a_typed_error(tmp_path: Path) -> None:
    """Given a malformed local info.json, When db list runs, Then provenance
    corruption is a typed DBINFO_MALFORMED error, not a silent local row."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fasta = tmp_path / "genes.fa"
    fasta.write_text(GENE_FASTA, encoding="utf-8")
    assert build("myamr", str(fasta), "--datadir", str(datadir)).exit_code == 0
    (datadir / "myamr" / "info.json").write_text("{not json", encoding="utf-8")

    result = runner.invoke(app, ["db", "list", "--datadir", str(datadir)])

    assert result.exit_code == 4
    envelope = json.loads(result.stderr.splitlines()[-1])
    assert envelope["code"] == "DBINFO_MALFORMED"
