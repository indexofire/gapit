"""`gapit db list` visibility for databases the registry does not know.

Databases installed into the datadir outside the provider catalog — `db
build` products (gene + cluster kinds) and manifest-less abricate-style
directories — must appear as ``local`` rows after the registry providers.
With no extras the output carries no trace of the extension (additive
proof). Real binaries through the CLI surface (test_cli_db_build.py
conventions).
"""

import json
import random
from pathlib import Path

from rich.console import Console
from typer.testing import CliRunner, Result

from gapit.cli import app
from gapit.db_ops import db_list_entries, db_list_table
from gapit.providers import REGISTRY

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
    product, When `db list`, Then both render as local rows after every
    registry provider — STATUS from the manifest record count, PROVIDER
    local, DBTYPE nucl, empty DESCRIPTION (no manifest note)."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fasta = tmp_path / "genes.fa"
    fasta.write_text(GENE_FASTA, encoding="utf-8")
    assert build("myamr", str(fasta), "--datadir", str(datadir)).exit_code == 0
    assert build("myloci", str(DATA / "bakta_style.gbk"), "--datadir", str(datadir)).exit_code == 0

    lines = list_rows(datadir)

    assert lines[0] == "NAME\tPROVIDER\tSTATUS\tDBTYPE\tDESCRIPTION"
    # Registry providers first (all available: nothing fetched), extras last
    # by name — deterministic ordering.
    assert lines[-2] == "myamr\tlocal\tinstalled (2)\tnucl\t"
    assert lines[-1] == "myloci\tlocal\tinstalled (2)\tnucl\t"
    assert all(line.split("\t")[1] != "local" for line in lines[1:-2])


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
    assert by_name["myloci"]["source"] == "local"
    assert by_name["myloci"]["kind"] == "cluster"
    registry_entries = [entry for entry in providers if "source" not in entry]
    assert registry_entries, "expected the real registry providers"
    assert all(entry["name"] in sorted(REGISTRY_NAMES) for entry in registry_entries)


def test_db_list_json_without_extras_carries_only_the_bundled_source_field(
    tmp_path: Path,
) -> None:
    """Given an empty datadir (no local extras), When `db list --json`, Then
    the only source fields in the document are the six bundled rows
    (alphabetical, source=bundled, not installed, vendor from bundled.json)
    — registry entries stay field-free and the extension remains additive."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = runner.invoke(app, ["db", "list", "--datadir", str(datadir), "--json"])

    assert result.exit_code == 0
    providers = json.loads(result.stdout)["providers"]
    with_source = [entry for entry in providers if "source" in entry]
    assert [(entry["name"], entry["source"]) for entry in with_source] == [
        ("ecoh", "bundled"),
        ("ecoli_dec", "bundled"),
        ("lm_doumith", "bundled"),
        ("ncbi", "bundled"),
        ("resfinder", "bundled"),
        ("upec_expec_vf", "bundled"),
    ]
    for entry in with_source:
        assert entry["installed"] is False
        assert "records" not in entry
        assert entry["dbtype"] == "nucl"
    assert with_source[1]["vendor"] == "gapit-curated (public-domain sources)"


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

    assert lines[-2] == "alegacy\tlocal\tinstalled (2)\tnucl\t"
    assert lines[-1] == "zproto\tlocal\tinstalled (1)\tprot\t"
