"""Bundled databases: discovery, materialization, and the DEC matrix.

The wheel ships five license-clean bundles — the typed ecoli_dec panel plus
four provider snapshots (ecoh, ncbi, resfinder, upec_expec_vf; public
domain / Apache-2.0 / BSD-3-Clause / MIT). These tests pin the
install-time-ready contract — discovery from the package data dir,
auto-materialization on first screen (idempotent, quiet-respecting,
fresh-datadir bootstrapping), the db list bundled → installed transition,
the per-snapshot screen smoke tests, the snapshot↔fetch differential
(resfinder, over a file:// provider URL), and the fetch-refresh path.
Real binaries through the CLI surface (test_cli_setupdb.py conventions).
"""

import json
import re
from dataclasses import replace
from pathlib import Path

import pytest
from typer.testing import CliRunner, Result

from gapit.bundled import bundled_databases, bundled_names, find_bundled
from gapit.cli import app
from gapit.dbcodec import decode_seqid
from gapit.fasta import iter_fasta
from gapit.providers import REGISTRY
from gapit.records import read_manifest

TYPING = Path(__file__).parent / "data" / "typing"
GOLDEN = Path(__file__).parent / "golden"

S1 = TYPING / "dec_s1_aggR_pic_uidA.fasta"
S2 = TYPING / "dec_s2_pic_astA_uidA.fasta"
S3 = TYPING / "dec_s3_stx2a_escV_aggR_uidA.fasta"

# The four provider snapshots: name -> committed record count (pins the
# wheel content against accidental regeneration).
SNAPSHOTS = {"ecoh": 597, "ncbi": 8373, "resfinder": 3206, "upec_expec_vf": 77}

# MUST-NOT guard: VFDB/CARD tags never ride the wheel (license audit
# 2026-10) — no VFG/VF id, no VFDB mention, no ARO accession, anywhere in
# any bundle's headers.
_RED_FLAG = re.compile(r"VF\d+|VFDB|VFG|ARO:", re.IGNORECASE)

RESFINDER_ZIP = Path(__file__).parent / "data" / "bundled_differential" / "resfinder" / "HEAD.zip"

runner = CliRunner()


def screen(sample: Path, datadir: Path, *extra: str) -> Result:
    return runner.invoke(
        app, ["screen", str(sample), "--db", "ecoli_dec", "--datadir", str(datadir), *extra]
    )


def phenotypes(sample: Path, datadir: Path, tmp_path: Path) -> dict[str, object]:
    """The DEC designation for one sample via the two-stage pipeline:
    screen to a table, then type it (gapit.typing_result/1)."""
    table = tmp_path / f"{sample.stem}.tsv"
    screened = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            "ecoli_dec",
            "--datadir",
            str(datadir),
            "-o",
            str(table),
            "--quiet",
        ],
    )
    assert screened.exit_code == 0, screened.stderr
    result = runner.invoke(
        app,
        ["typing", str(table), "--datadir", str(datadir), "--format", "json", "--quiet"],
    )
    assert result.exit_code == 0, result.stderr
    document = json.loads(result.stdout)
    assert document["db"] == "ecoli_dec"
    return {name: call["phenotype"] for name, call in document["files"][0]["phenotypes"].items()}


def sample_from_bundle(name: str, destination: Path, count: int = 3) -> set[str]:
    """Write the first ``count`` bundled records as one contig each; return
    their decoded gene names (the expected screen calls)."""
    bundled = find_bundled(name)
    assert bundled is not None
    genes: set[str] = set()
    with destination.open("w", encoding="utf-8") as handle:
        for record in [
            r for _, r in zip(range(count), iter_fasta(bundled.sequences_path), strict=False)
        ]:
            handle.write(f">{record.id} sampled\n{record.sequence}\n")
            genes.add(decode_seqid(record.id, default_db=name).gene)
    return genes


def screen_db(sample: Path, datadir: Path, name: str, *extra: str) -> Result:
    return runner.invoke(
        app, ["screen", str(sample), "--db", name, "--datadir", str(datadir), *extra]
    )


def test_bundled_discovery_finds_all_five() -> None:
    """Given the wheel's data dir, When bundled databases are discovered,
    Then exactly the five license-clean bundles ship in name order, each
    metadata parses (vendor/dbtype/snapshot date), ecoli_dec alone carries
    the dual-scheme typing document, and the provider snapshots do not."""
    assert bundled_names() == ["ecoh", "ecoli_dec", "ncbi", "resfinder", "upec_expec_vf"]
    by_name = {database.name: database for database in bundled_databases()}
    ecoli_dec = by_name["ecoli_dec"]
    assert ecoli_dec.metadata.vendor == "gapit-curated (public-domain sources)"
    assert ecoli_dec.metadata.dbtype == "nucl"
    assert ecoli_dec.metadata.snapshotted == "2026-10-02"
    assert ecoli_dec.typing_path is not None
    vendors = {
        "ecoh": "Holt lab (srst2)",
        "ncbi": "NCBI",
        "resfinder": "DTU CGE",
        "upec_expec_vf": "FordeGenomics",
    }
    for name, vendor in vendors.items():
        assert by_name[name].metadata.vendor == vendor
        assert by_name[name].metadata.dbtype == "nucl"
        assert by_name[name].typing_path is None
    assert find_bundled("ecoli_dec") is not None
    assert find_bundled("ncbi") is not None
    assert find_bundled("card") is None


def test_bundled_headers_carry_no_red_flag_tags() -> None:
    """Given every bundled panel, When every header is inspected, Then no VF
    tag, VFDB marker, or ARO accession appears anywhere — the DEC lesson
    from the license audit (2 VFDB records hid past header-spotting) is
    locked mechanically across all five bundles."""
    for bundled in bundled_databases():
        text = bundled.sequences_path.read_text(encoding="utf-8")
        headers = [line for line in text.splitlines() if line.startswith(">")]
        if bundled.name == "ecoli_dec":
            assert len(headers) == 17
        else:
            assert len(headers) == SNAPSHOTS[bundled.name]
        for header in headers:
            assert _RED_FLAG.search(header) is None, header


def test_screen_auto_materializes_on_fresh_empty_datadir(tmp_path: Path) -> None:
    """Given a fresh EMPTY datadir and no prior fetch, When screening
    against the bundled ecoli_dec, Then the database materializes into the
    datadir first (one stderr note), the screen finds the panel genes, and
    the manifest certifies a bundled-source typed gene database."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = screen(S1, datadir)

    assert result.exit_code == 0, result.stderr
    assert f"gapit: materializing bundled database ecoli_dec (17 records) into {datadir}" in (
        result.stderr
    )
    db_dir = datadir / "ecoli_dec"
    manifest = read_manifest(db_dir / "gapit-manifest.json")
    assert manifest.n_records == 17
    assert manifest.source == "bundled"
    assert manifest.typing_schema == "gapit.typing/2"
    assert (db_dir / "typing.json").is_file()
    assert (db_dir / "sequences.nin").is_file()
    genes = {line.split("\t")[5] for line in result.stdout.splitlines()[1:]}
    assert genes == {"aggR", "pic", "uidA"}


def test_screen_is_idempotent_and_quiet_suppresses_the_note(tmp_path: Path) -> None:
    """Given a materialized ecoli_dec, When screening again, Then no rebuild
    happens (no materializing note, manifest bytes unchanged); with --quiet
    even the first materialization note is silent."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    first = screen(S1, datadir)
    assert first.exit_code == 0, first.stderr
    manifest_path = datadir / "ecoli_dec" / "gapit-manifest.json"
    before = manifest_path.read_bytes()

    second = screen(S2, datadir)

    assert second.exit_code == 0, second.stderr
    assert "materializing" not in second.stderr
    assert manifest_path.read_bytes() == before

    other = tmp_path / "quiet-dd"
    other.mkdir()
    silent = screen(S1, other, "--quiet")
    assert silent.exit_code == 0, silent.stderr
    assert silent.stderr == ""
    assert (other / "ecoli_dec" / "gapit-manifest.json").is_file()


def test_screen_bootstraps_a_missing_datadir_for_bundled(tmp_path: Path) -> None:
    """Given a --datadir that does not exist at all (fresh machine), When
    screening against a bundled database, Then the datadir is created and
    the database materialized — zero-network first use; a NON-bundled db
    keeps the typed DATADIR_NOT_FOUND error."""
    missing = tmp_path / "does" / "not" / "exist"
    result = screen(S1, missing)
    assert result.exit_code == 0, result.stderr
    assert (missing / "ecoli_dec" / "gapit-manifest.json").is_file()

    other = tmp_path / "also-missing"
    failure = runner.invoke(app, ["screen", str(S1), "--db", "card", "--datadir", str(other)])
    assert failure.exit_code == 4
    assert json.loads(failure.stderr)["code"] == "DATADIR_NOT_FOUND"


def test_db_list_bundled_to_installed_transition(tmp_path: Path) -> None:
    """Given an empty datadir, When `db list`, Then ecoli_dec renders as a
    bundled row (vendor from bundled.json, STATUS bundled) between the
    registry rows and any local extras; after one screen the same row reads
    installed (17) and the JSON entry flips installed/source stays bundled."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    before = runner.invoke(app, ["db", "list", "--datadir", str(datadir)])
    assert before.exit_code == 0
    lines = before.stdout.splitlines()
    bundled_rows = [line for line in lines if line.startswith("ecoli_dec\t")]
    assert bundled_rows == [
        "ecoli_dec\tgapit-curated (public-domain sources)\tbundled\tnucl\t"
        "Diarrheagenic E. coli marker panel (GB 4789.6 + risk-monitoring designation)"
    ]
    bundled_five = set(bundled_names())
    first_bundled = next(
        index
        for index, line in enumerate(lines[1:], start=1)
        if line.split("\t")[0] in bundled_five
    )
    registry_rows = lines[1:first_bundled]
    assert registry_rows and all(
        row.split("\t")[2] == "available" for row in registry_rows if row.split("\t")[2] != "local"
    )

    assert screen(S1, datadir).exit_code == 0
    after = runner.invoke(app, ["db", "list", "--datadir", str(datadir)])
    assert after.exit_code == 0
    assert (
        "ecoli_dec\tgapit-curated (public-domain sources)\tinstalled (17)\tnucl\t"
        "Diarrheagenic E. coli marker panel (GB 4789.6 + risk-monitoring designation)"
        in after.stdout
    )

    as_json = runner.invoke(app, ["db", "list", "--datadir", str(datadir), "--json"])
    entry = {
        e["name"]: e for e in json.loads(as_json.stdout)["providers"] if e["name"] == "ecoli_dec"
    }["ecoli_dec"]
    assert entry["source"] == "bundled"
    assert entry["installed"] is True
    assert entry["records"] == 17


def test_dec_designation_matrix_against_the_bundled_database(tmp_path: Path) -> None:
    """Given the bundled ecoli_dec (typed gapit.typing/2), When the three
    DEC matrix samples run the two-stage pipeline (screen -o table, then
    typing table), Then the dual-scheme calls hold — the equivalence lock
    with the inline engine's former calls: aggR+pic → EAEC/EAEC, pic+astA
    → EAEC/non-DEC (the headline scheme divergence), stx2a+escV hybrid →
    EHEC/EHEC."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    assert screen(S1, datadir).exit_code == 0

    assert phenotypes(S1, datadir, tmp_path) == {"gb4789_6": "EAEC", "risk_monitoring": "EAEC"}
    assert phenotypes(S2, datadir, tmp_path) == {"gb4789_6": "EAEC", "risk_monitoring": "non-DEC"}
    assert phenotypes(S3, datadir, tmp_path) == {"gb4789_6": "EHEC", "risk_monitoring": "EHEC"}


def test_typed_screen_golden_hybrid_and_headline(tmp_path: Path) -> None:
    """Given the bundled ecoli_dec and the hybrid + headline samples, When
    screened as the DEFAULT format (multi-file, --nopath), Then the output
    is the pure frozen 15-column abricate table — typed and untyped gene
    databases screen byte-identically (the inverted former
    phenotype-column golden; designation now renders via `gapit typing`)."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    result = runner.invoke(
        app,
        [
            "screen",
            str(S3),
            str(S2),
            "--db",
            "ecoli_dec",
            "--datadir",
            str(datadir),
            "--nopath",
            "--quiet",
        ],
    )
    assert result.exit_code == 0, result.stderr
    assert result.stdout == (GOLDEN / "ecoli_dec_typed.tsv").read_text(encoding="utf-8")
    for line in result.stdout.splitlines():
        assert len(line.split("\t")) == 15


def test_typing_golden_hybrid_and_headline(tmp_path: Path) -> None:
    """Given the screened hybrid + headline tables, When typed as the
    DEFAULT TSV, Then the seven-column designation output is byte-identical
    to the committed golden — the hybrid's EHEC/EHEC beside the headline
    EAEC/non-DEC divergence."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    table = tmp_path / "screen.tsv"
    screened = runner.invoke(
        app,
        [
            "screen",
            str(S3),
            str(S2),
            "--db",
            "ecoli_dec",
            "--datadir",
            str(datadir),
            "--nopath",
            "-o",
            str(table),
            "--quiet",
        ],
    )
    assert screened.exit_code == 0, screened.stderr
    result = runner.invoke(app, ["typing", str(table), "--datadir", str(datadir), "--quiet"])
    assert result.exit_code == 0, result.stderr
    assert result.stdout == (GOLDEN / "ecoli_dec_typing.tsv").read_text(encoding="utf-8")


# ------------------------------------------------- provider snapshot bundles --


@pytest.mark.parametrize(("db_name", "db_records"), sorted(SNAPSHOTS.items()))
def test_snapshot_materializes_and_screens_zero_network(
    tmp_path: Path, db_name: str, db_records: int
) -> None:
    """Given a fresh EMPTY datadir and a sample derived from the bundle's own
    sequences, When screening --db <snapshot>, Then the bundle materializes
    (note carries the committed record count), the manifest certifies
    source=bundled with no typing document, and every sampled gene is called
    at 100% identity and coverage — the zero-network headline proof."""
    sample = tmp_path / "sample.fa"
    expected = sample_from_bundle(db_name, sample)
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = screen_db(sample, datadir, db_name)

    assert result.exit_code == 0, result.stderr
    assert f"materializing bundled database {db_name} ({db_records} records)" in result.stderr
    manifest = read_manifest(datadir / db_name / "gapit-manifest.json")
    assert manifest.n_records == db_records
    assert manifest.source == "bundled"
    assert manifest.typing_schema is None
    assert not (datadir / db_name / "typing.json").exists()
    rows = [line.split("\t") for line in result.stdout.splitlines()[1:]]
    called = {row[5] for row in rows}
    assert expected <= called
    assert all(row[9] == "100.00" for row in rows if row[5] in expected)
    assert all(row[10] == "100.00" for row in rows if row[5] in expected)


def test_db_list_merges_registry_bundled_names_into_one_row(tmp_path: Path) -> None:
    """Given an empty datadir and the REAL registry, When `db list`, Then
    each snapshot name renders exactly ONCE — as a bundled row (STATUS
    bundled, vendor from bundled.json) inside the five-row alphabetical
    bundled section, never as a duplicate registry available row — while a
    NON-bundled registry row keeps its exact registry bytes; after one
    upec_expec_vf screen that row flips to installed (77) and the JSON
    entry carries source=bundled."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    before = runner.invoke(app, ["db", "list", "--datadir", str(datadir)])
    assert before.exit_code == 0
    lines = before.stdout.splitlines()
    assert lines[-5:] == [
        "ecoh\tHolt lab (srst2)\tbundled\tnucl\tE. coli O and H antigens (srst2 EcOH)",
        "ecoli_dec\tgapit-curated (public-domain sources)\tbundled\tnucl\t"
        "Diarrheagenic E. coli marker panel (GB 4789.6 + risk-monitoring designation)",
        "ncbi\tNCBI\tbundled\tnucl\tNCBI AMRFinderPlus (reference finder) curated AMR",
        "resfinder\tDTU CGE\tbundled\tnucl\tCGE ResFinder acquired resistance genes",
        "upec_expec_vf\tFordeGenomics\tbundled\tnucl\tUPEC/ExPEC virulence genes (FordeGenomics)",
    ]
    for name in SNAPSHOTS:
        assert sum(line.startswith(f"{name}\t") for line in lines) == 1, name
    assert (
        "card\tMcMaster University\tavailable\tnucl\tCARD protein homolog resistance models"
        in (lines[1:-5])
    )

    sample = tmp_path / "sample.fa"
    sample_from_bundle("upec_expec_vf", sample)
    assert screen_db(sample, datadir, "upec_expec_vf", "--quiet").exit_code == 0

    after = runner.invoke(app, ["db", "list", "--datadir", str(datadir)])
    assert after.exit_code == 0
    assert (
        "upec_expec_vf\tFordeGenomics\tinstalled (77)\tnucl\t"
        "UPEC/ExPEC virulence genes (FordeGenomics)" in after.stdout
    )

    as_json = runner.invoke(app, ["db", "list", "--datadir", str(datadir), "--json"])
    providers = json.loads(as_json.stdout)["providers"]
    by_name = [entry["name"] for entry in providers]
    assert by_name.count("upec_expec_vf") == 1
    entry = {entry["name"]: entry for entry in providers}["upec_expec_vf"]
    assert entry["source"] == "bundled"
    assert entry["installed"] is True
    assert entry["records"] == 77


def test_resfinder_snapshot_is_byte_identical_to_a_fetch_build(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the wheel's resfinder snapshot and a file:// copy of the exact
    upstream HEAD.zip it was produced from, When one datadir is built by
    `db fetch` (the provider pipeline) and another by bundled
    materialization, Then both sequences files are BYTE-identical and carry
    the same record count — the snapshot is precisely what fetch installs."""
    provider = replace(REGISTRY["resfinder"], source_urls=(RESFINDER_ZIP.as_uri(),))
    monkeypatch.setattr("gapit.db_ops.REGISTRY", {"resfinder": provider})

    fetch_dd = tmp_path / "fetch-dd"
    fetch_dd.mkdir()
    fetched = runner.invoke(
        app, ["db", "fetch", "resfinder", "--datadir", str(fetch_dd), "--quiet"]
    )
    assert fetched.exit_code == 0, fetched.stderr
    assert json.loads(fetched.stdout)["records"] == SNAPSHOTS["resfinder"]

    sample = tmp_path / "sample.fa"
    sample_from_bundle("resfinder", sample)
    bundle_dd = tmp_path / "bundle-dd"
    bundle_dd.mkdir()
    assert screen_db(sample, bundle_dd, "resfinder", "--quiet").exit_code == 0

    fetched_sequences = fetch_dd / "resfinder" / "sequences"
    bundled_sequences = bundle_dd / "resfinder" / "sequences"
    assert bundled_sequences.read_bytes() == fetched_sequences.read_bytes()


def test_db_fetch_remains_the_refresh_path_for_a_bundled_name(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a materialized resfinder bundle, When `db fetch resfinder`
    (file:// provider), Then the first call refuses (DB_ALREADY_EXISTS — a
    materialized bundle is an installed database) and --force overwrites it
    with the fresh upstream build: the manifest drops its bundled source
    stamp while `db list` keeps the row in the bundled section (wheel
    membership) as installed — the snapshot is point-in-time, fetch is the
    freshness path."""
    provider = replace(REGISTRY["resfinder"], source_urls=(RESFINDER_ZIP.as_uri(),))
    monkeypatch.setattr("gapit.db_ops.REGISTRY", {"resfinder": provider})
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    sample = tmp_path / "sample.fa"
    sample_from_bundle("resfinder", sample)
    assert screen_db(sample, datadir, "resfinder", "--quiet").exit_code == 0
    manifest_path = datadir / "resfinder" / "gapit-manifest.json"
    assert read_manifest(manifest_path).source == "bundled"

    refusal = runner.invoke(app, ["db", "fetch", "resfinder", "--datadir", str(datadir)])
    assert refusal.exit_code == 4
    assert json.loads(refusal.stderr)["code"] == "DB_ALREADY_EXISTS"

    refresh = runner.invoke(
        app, ["db", "fetch", "resfinder", "--datadir", str(datadir), "--force", "--quiet"]
    )
    assert refresh.exit_code == 0, refresh.stderr
    assert read_manifest(manifest_path).source is None

    listing = runner.invoke(app, ["db", "list", "--datadir", str(datadir), "--json"])
    entry = {entry["name"]: entry for entry in json.loads(listing.stdout)["providers"]}["resfinder"]
    assert entry["source"] == "bundled"
    assert entry["installed"] is True
    assert entry["records"] == SNAPSHOTS["resfinder"]
