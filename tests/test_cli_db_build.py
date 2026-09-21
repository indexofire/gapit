"""CLI integration tests: `gapit db build NAME FASTA` (custom databases).

Real binaries (makeblastdb/blastn/minimap2 from the pixi env), synthetic
fixtures only, everything local to tmp_path — the test_cli_db_fetch.py and
test_gapit_db_e2e.py conventions. Header auto-detection (plain / abricate
``~~~`` / ``gapit|``), the --tsv metadata merge, dbtype resolution, the
DB_ALREADY_EXISTS/--force gate, and one full E2E screen against a built
database per header kind.
"""

import json
import random
from pathlib import Path

from pydantic import TypeAdapter
from typer.testing import CliRunner, Result

from gapit.cli import app
from gapit.errors import ErrorEnvelope
from gapit.records import Record, read_manifest, read_records

DB = "myamr"
OTHER_DB = "myamr2"

# Distinct seeds share no 11-mer in practice (test_gapit_db_e2e.py note), so
# the synthetic genes cannot cross-match under blastn.
SEQ_A = "".join(random.Random(1).choice("ACGT") for _ in range(240))  # demov2
SEQ_B = "".join(random.Random(2).choice("ACGT") for _ in range(240))  # demov3
SEQ_C = "".join(random.Random(3).choice("ACGT") for _ in range(240))  # tilde gene
SEQ_D = "".join(random.Random(4).choice("ACGT") for _ in range(240))  # gapit gene
AMINO = "".join(
    random.Random(5).choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(80)
)  # protein-looking

PLAIN_FASTA = f">demov2 demo beta-lactamase variant 2\n{SEQ_A}\n>demov3 demo efflux pump\n{SEQ_B}\n"
TILDE_SEQID = "olddb~~~tilde1~~~SYN-001~~~SULFONAMIDE"
TILDE_FASTA = f">{TILDE_SEQID} tilde demo protein\n{SEQ_C}\n"
GAPIT_SEQID = "gapit|db=olddb|gene=tagged1|acc=SYN-002|func=ampicillin;gentamicin"
GAPIT_FASTA = f">{GAPIT_SEQID} tagged demo protein\n{SEQ_D}\n"

runner = CliRunner()
envelope_adapter = TypeAdapter(ErrorEnvelope)


def last_envelope(stderr: str) -> ErrorEnvelope:
    """Parse the last non-empty stderr line as the gapit.error/1 envelope."""
    lines = [line for line in stderr.splitlines() if line.strip()]
    assert lines, "expected an error envelope on stderr"
    return envelope_adapter.validate_json(lines[-1])


def build(*extra: str) -> Result:
    return runner.invoke(app, ["db", "build", *extra])


def records_of(datadir: Path, name: str = DB) -> list[Record]:
    return list(read_records(datadir / name / "records.jsonl"))


def screen_rows(datadir: Path, query: Path, name: str = DB) -> list[list[str]]:
    """Run a real screen and return the TSV data rows (header stripped)."""
    result = runner.invoke(app, ["screen", str(query), "--db", name, "--datadir", str(datadir)])
    assert result.exit_code == 0, result.stderr
    return [line.split("\t") for line in result.stdout.splitlines() if not line.startswith("#")]


def test_db_build_plain_fasta_creates_database_and_screens(tmp_path: Path) -> None:
    """Given a plain two-gene FASTA, When `db build NAME FASTA`, Then exit 0
    with a fetch-shaped JSON receipt, every native artifact exists (sequences,
    .nin, manifest, records.jsonl) and no .mmi is built, the manifest
    certifies a local build, and screening a 200/240 bp substring of one gene
    returns exactly that hit with gene and product decoded."""
    fasta = tmp_path / "my_genes.fa"
    fasta.write_text(PLAIN_FASTA, encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build(DB, str(fasta), "--datadir", str(datadir))

    assert result.exit_code == 0, result.stderr
    db_dir = datadir / DB
    assert json.loads(result.stdout) == {
        "db": DB,
        "records": 2,
        "dbtype": "nucl",
        "destination": str(db_dir),
    }
    for artifact in ("sequences", "sequences.nin", "gapit-manifest.json"):
        assert (db_dir / artifact).is_file(), artifact
    assert not (db_dir / "sequences.mmi").exists()
    manifest = read_manifest(db_dir / "gapit-manifest.json")
    assert manifest.source_urls == ("local",)
    assert manifest.upstream_version == ""
    assert manifest.dbtype == "nucl"

    query = tmp_path / "query.fa"
    query.write_text(f">contig1\n{SEQ_A[:200]}\n", encoding="utf-8")
    (row,) = screen_rows(datadir, query)
    assert row[5] == "demov2"  # GENE
    assert row[13] == "demo beta-lactamase variant 2"  # PRODUCT
    assert row[12] == ""  # ACCESSION: plain header carries none


def test_db_build_abricate_tilde_headers(tmp_path: Path) -> None:
    """Given an abricate `~~~` FASTA, When built and screened with the full
    gene as query, Then the TSV row shows the gene/accession/function fields
    extracted from the header and the product from the description, with the
    Record db field retargeted to NAME."""
    fasta = tmp_path / "tilde.fa"
    fasta.write_text(TILDE_FASTA, encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build(DB, str(fasta), "--datadir", str(datadir))

    assert result.exit_code == 0, result.stderr
    (record,) = records_of(datadir)
    assert record.db == DB
    assert record.gene == "tilde1"
    assert record.accession == "SYN-001"
    assert record.function == ("SULFONAMIDE",)
    assert record.product == "tilde demo protein"
    assert record.source_id == TILDE_SEQID

    query = tmp_path / "query.fa"
    query.write_text(f">contig1\n{SEQ_C}\n", encoding="utf-8")
    (row,) = screen_rows(datadir, query)
    assert row[5] == "tilde1"
    assert row[12] == "SYN-001"
    assert row[13] == "tilde demo protein"
    assert row[14] == "SULFONAMIDE"


def test_db_build_gapit_headers_roundtrip(tmp_path: Path) -> None:
    """Given a gapit| tagged FASTA, When built and screened, Then gene,
    accession and the multi-class function decode back through screening
    (roundtrip) and source_id keeps the original id token."""
    fasta = tmp_path / "tagged.fa"
    fasta.write_text(GAPIT_FASTA, encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build(DB, str(fasta), "--datadir", str(datadir))

    assert result.exit_code == 0, result.stderr
    (record,) = records_of(datadir)
    assert record.db == DB
    assert record.gene == "tagged1"
    assert record.accession == "SYN-002"
    assert record.function == ("ampicillin", "gentamicin")
    assert record.source_id.startswith("gapit|db=olddb|")

    query = tmp_path / "query.fa"
    query.write_text(f">contig1\n{SEQ_D}\n", encoding="utf-8")
    (row,) = screen_rows(datadir, query)
    assert row[5] == "tagged1"
    assert row[12] == "SYN-002"
    assert row[14] == "ampicillin;gentamicin"


def test_db_build_mixed_header_file(tmp_path: Path) -> None:
    """Given one FASTA mixing plain, `~~~` and gapit| headers, When built,
    Then each record is parsed per its own header kind (db = NAME for all)."""
    fasta = tmp_path / "mixed.fa"
    fasta.write_text(f"{PLAIN_FASTA}{TILDE_FASTA}{GAPIT_FASTA}", encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build(DB, str(fasta), "--datadir", str(datadir))

    assert result.exit_code == 0, result.stderr
    parsed = records_of(datadir)
    assert [record.gene for record in parsed] == ["demov2", "demov3", "tilde1", "tagged1"]
    assert all(record.db == DB for record in parsed)
    assert parsed[0].accession == ""
    assert parsed[2].accession == "SYN-001"
    assert parsed[3].function == ("ampicillin", "gentamicin")


def test_db_build_description_flag_fills_product(tmp_path: Path) -> None:
    """Given a plain FASTA record with no description text, When built with
    --description TEXT, Then that text becomes the product; without the flag
    the product falls back to the gene name."""
    fasta = tmp_path / "bare.fa"
    fasta.write_text(f">baregene\n{SEQ_A}\n", encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    with_flag = build("flagged", str(fasta), "--datadir", str(datadir), "--description", "custom")
    assert with_flag.exit_code == 0, with_flag.stderr
    (record,) = records_of(datadir, "flagged")
    assert record.product == "custom"

    without = build("unflagged", str(fasta), "--datadir", str(datadir))
    assert without.exit_code == 0, without.stderr
    (record,) = records_of(datadir, "unflagged")
    assert record.product == "baregene"


def test_db_build_tsv_merges_accession_and_function(tmp_path: Path) -> None:
    """Given a --tsv with accession+function columns, When built, Then both
    fields overwrite the FASTA-parsed values for matching genes and the
    semicolon-separated classes screen as one RESISTANCE string."""
    fasta = tmp_path / "my_genes.fa"
    fasta.write_text(PLAIN_FASTA, encoding="utf-8")
    tsv = tmp_path / "meta.tsv"
    tsv.write_text(
        "gene\taccession\tfunction\tnote\ndemov2\tSYN-X\tampicillin;gentamicin\tignored\n",
        encoding="utf-8",
    )
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build(DB, str(fasta), "--datadir", str(datadir), "--tsv", str(tsv))

    assert result.exit_code == 0, result.stderr
    by_gene = {record.gene: record for record in records_of(datadir)}
    assert by_gene["demov2"].accession == "SYN-X"
    assert by_gene["demov2"].function == ("ampicillin", "gentamicin")
    assert by_gene["demov3"].accession == ""  # no TSV row: untouched

    query = tmp_path / "query.fa"
    query.write_text(f">contig1\n{SEQ_A}\n", encoding="utf-8")
    (row,) = screen_rows(datadir, query)
    assert row[12] == "SYN-X"
    assert row[14] == "ampicillin;gentamicin"


def test_db_build_tsv_partial_columns_merge_only_present(tmp_path: Path) -> None:
    """Given a --tsv whose header carries gene+function but no accession
    column, When built, Then only the function is overwritten."""
    fasta = tmp_path / "my_genes.fa"
    fasta.write_text(PLAIN_FASTA, encoding="utf-8")
    tsv = tmp_path / "meta.tsv"
    tsv.write_text("gene\tfunction\ndemov2\tvirulence\n", encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build(DB, str(fasta), "--datadir", str(datadir), "--tsv", str(tsv))

    assert result.exit_code == 0, result.stderr
    (record,) = (r for r in records_of(datadir) if r.gene == "demov2")
    assert record.function == ("virulence",)
    assert record.accession == ""


def test_db_build_tsv_duplicate_gene_first_row_wins(tmp_path: Path) -> None:
    """Given a --tsv with two rows for one gene, When built, Then the FIRST
    row wins and a warning lands on stderr (--quiet silences it)."""
    fasta = tmp_path / "my_genes.fa"
    fasta.write_text(f">demov2 demo\n{SEQ_A}\n", encoding="utf-8")
    tsv = tmp_path / "dups.tsv"
    tsv.write_text(
        "gene\taccession\tfunction\ndemov2\tFIRST\ttetracycline\ndemov2\tSECOND\tvirulence\n",
        encoding="utf-8",
    )
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    loud = build(DB, str(fasta), "--datadir", str(datadir), "--tsv", str(tsv))
    assert loud.exit_code == 0, loud.stderr
    (record,) = records_of(datadir)
    assert record.accession == "FIRST"
    assert record.function == ("tetracycline",)
    assert "demov2" in loud.stderr
    assert "WARNING" in loud.stderr

    quiet = build(OTHER_DB, str(fasta), "--datadir", str(datadir), "--tsv", str(tsv), "--quiet")
    assert quiet.exit_code == 0, quiet.stderr
    assert quiet.stderr == ""


def test_db_build_tsv_gene_absent_from_fasta_warns_and_skips(tmp_path: Path) -> None:
    """Given a --tsv row naming a gene the FASTA does not carry, When built,
    Then the build succeeds, the row is skipped, and a warning names it."""
    fasta = tmp_path / "my_genes.fa"
    fasta.write_text(f">demov2 demo\n{SEQ_A}\n", encoding="utf-8")
    tsv = tmp_path / "ghost.tsv"
    tsv.write_text("gene\taccession\tfunction\nghost\tSYN-G\tvirulence\n", encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build(DB, str(fasta), "--datadir", str(datadir), "--tsv", str(tsv))

    assert result.exit_code == 0, result.stderr
    (record,) = records_of(datadir)
    assert record.gene == "demov2"
    assert "ghost" in result.stderr


def test_db_build_tsv_missing_gene_column_is_malformed(tmp_path: Path) -> None:
    """Given a --tsv whose header lacks the mandatory gene column, When built,
    Then exit 5 with InputError METADATA_MALFORMED and {file, needed} context."""
    fasta = tmp_path / "my_genes.fa"
    fasta.write_text(PLAIN_FASTA, encoding="utf-8")
    tsv = tmp_path / "bad.tsv"
    tsv.write_text("accession\tfunction\nSYN-X\tvirulence\n", encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build(DB, str(fasta), "--datadir", str(datadir), "--tsv", str(tsv))

    assert result.exit_code == 5
    envelope = last_envelope(result.stderr)
    assert envelope.code == "METADATA_MALFORMED"
    assert envelope.context["file"] == str(tsv)
    assert envelope.context["needed"] == "gene"


def test_db_build_missing_fasta_is_input_not_found(tmp_path: Path) -> None:
    """Given a FASTA path that does not exist, When built, Then exit 5 with
    INPUT_NOT_FOUND naming the file."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    result = build(DB, str(tmp_path / "nope.fa"), "--datadir", str(datadir))
    assert result.exit_code == 5
    envelope = last_envelope(result.stderr)
    assert envelope.code == "INPUT_NOT_FOUND"
    assert envelope.context["file"] == str(tmp_path / "nope.fa")


def test_db_build_missing_tsv_is_input_not_found(tmp_path: Path) -> None:
    """Given a --tsv path that does not exist (but a valid FASTA), When built,
    Then exit 5 with INPUT_NOT_FOUND naming the TSV."""
    fasta = tmp_path / "my_genes.fa"
    fasta.write_text(PLAIN_FASTA, encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    result = build(DB, str(fasta), "--datadir", str(datadir), "--tsv", str(tmp_path / "nope.tsv"))
    assert result.exit_code == 5
    envelope = last_envelope(result.stderr)
    assert envelope.code == "INPUT_NOT_FOUND"
    assert envelope.context["file"] == str(tmp_path / "nope.tsv")


def test_db_build_invalid_fasta_propagates_typed(tmp_path: Path) -> None:
    """Given a structurally invalid FASTA (content before the first header),
    When built, Then the reader's INVALID_FASTA envelope propagates (exit 5)."""
    fasta = tmp_path / "bad.fa"
    fasta.write_text("not a fasta\n", encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    result = build(DB, str(fasta), "--datadir", str(datadir))
    assert result.exit_code == 5
    envelope = last_envelope(result.stderr)
    assert envelope.code == "INVALID_FASTA"


def test_db_build_prot_dbtype_builds_pin(tmp_path: Path) -> None:
    """Given a protein FASTA and --dbtype prot, When built, Then the receipt
    declares prot and the BLAST index is .pin."""
    fasta = tmp_path / "prot.faa"
    fasta.write_text(f">prot1 synthetic protein\n{AMINO}\n", encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build(DB, str(fasta), "--datadir", str(datadir), "--dbtype", "prot")

    assert result.exit_code == 0, result.stderr
    assert json.loads(result.stdout)["dbtype"] == "prot"
    assert (datadir / DB / "sequences.pin").is_file()


def test_db_build_auto_detects_protein_without_flag(tmp_path: Path) -> None:
    """Given a protein-looking FASTA with NO --dbtype, When built, Then the
    mol_type heuristic picks prot (receipt + manifest agree)."""
    fasta = tmp_path / "prot.faa"
    fasta.write_text(f">prot1 synthetic protein\n{AMINO}\n", encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build(DB, str(fasta), "--datadir", str(datadir))

    assert result.exit_code == 0, result.stderr
    assert json.loads(result.stdout)["dbtype"] == "prot"
    assert read_manifest(datadir / DB / "gapit-manifest.json").dbtype == "prot"


def test_db_build_refuses_existing_database_then_force_rebuilds(tmp_path: Path) -> None:
    """Given an already-built database, When built again without --force,
    Then exit 4 with DB_ALREADY_EXISTS; When built with --force, Then exit 0
    with a fresh receipt over the same destination."""
    fasta = tmp_path / "my_genes.fa"
    fasta.write_text(PLAIN_FASTA, encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    assert build(DB, str(fasta), "--datadir", str(datadir)).exit_code == 0

    refused = build(DB, str(fasta), "--datadir", str(datadir))
    assert refused.exit_code == 4
    assert refused.stdout == ""
    envelope = last_envelope(refused.stderr)
    assert envelope.code == "DB_ALREADY_EXISTS"
    assert envelope.context["db"] == DB

    forced = build(DB, str(fasta), "--datadir", str(datadir), "--force")
    assert forced.exit_code == 0, forced.stderr
    assert json.loads(forced.stdout)["records"] == 2


def test_db_build_rejects_names_that_escape_the_datadir(tmp_path: Path) -> None:
    """Given names that would redirect the destination outside the datadir
    (absolute, ../ traversal, separator, dot), When built, Then exit 2 with a
    USAGE_ERROR envelope carrying the name in context and nothing is written
    anywhere (datadir stays empty, no ../escape sibling)."""
    fasta = tmp_path / "my_genes.fa"
    fasta.write_text(PLAIN_FASTA, encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    for name in ("/tmp/evil", "../escape", "a/b", "."):
        result = build(name, str(fasta), "--datadir", str(datadir))
        assert result.exit_code == 2, name
        envelope = last_envelope(result.stderr)
        assert envelope.code == "USAGE_ERROR"
        assert envelope.context["name"] == name
        assert envelope.message == (
            f"database name must be a plain name without path separators: {name!r}"
        )
    assert list(datadir.iterdir()) == []
    assert not (tmp_path / "escape").exists()


def test_db_build_allows_dash_underscore_dot_names(tmp_path: Path) -> None:
    """Given a plain name carrying dashes, underscores and an inside dot,
    When built, Then exit 0 — anything without separators stays allowed."""
    fasta = tmp_path / "my_genes.fa"
    fasta.write_text(f">demov2 demo\n{SEQ_A}\n", encoding="utf-8")
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = build("my-db_1.x", str(fasta), "--datadir", str(datadir))

    assert result.exit_code == 0, result.stderr
    assert (datadir / "my-db_1.x" / "gapit-manifest.json").is_file()


def test_db_build_creates_missing_datadir(tmp_path: Path) -> None:
    """Given a --datadir whose nested path does not exist yet, When built,
    Then the datadir is created (mirroring `db fetch` bootstrap) and the
    database lands under it."""
    fasta = tmp_path / "my_genes.fa"
    fasta.write_text(PLAIN_FASTA, encoding="utf-8")
    nested = tmp_path / "fresh" / "machine" / "db"

    result = build(DB, str(fasta), "--datadir", str(nested))

    assert result.exit_code == 0, result.stderr
    assert json.loads(result.stdout)["destination"] == str(nested / DB)
    assert (nested / DB / "gapit-manifest.json").is_file()
