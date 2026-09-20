"""Integration tests: the real any2fasta -> blast pipeline over committed contig
fixtures against the tinyamr database (offline, uses the pixi env's binaries)."""

import shutil
from pathlib import Path

import pytest

from gapit.blast import screen_file
from gapit.db import Database, discover_databases, make_blast_db
from gapit.errors import InputError
from gapit.report import ScreeningParams

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"
CONTIGS = Path(__file__).parent / "data" / "contigs"
PARAMS = ScreeningParams(db="tinyamr")


@pytest.fixture()
def tinyamr(tmp_path: Path) -> Database:
    """Fresh datadir with a built tinyamr index, one per test."""
    datadir = tmp_path / "datadir"
    shutil.copytree(FIXTURE_DB_DIR, datadir)
    (database,) = discover_databases(datadir)
    make_blast_db(database.sequences_path, database.name)
    return database


def test_full_length_exact_match(tinyamr: Database) -> None:
    """Given a contig that IS the tetA gene, When screened, Then one 100%/100%
    hit with an all-'=' map and the cleaned product."""
    report = screen_file(CONTIGS / "full.fa", tinyamr, PARAMS, dbtype="nucl")
    assert report.file == str(CONTIGS / "full.fa")
    (hit,) = report.hits
    assert hit.sequence == "contig1"
    assert (hit.start, hit.end) == (1, 79)
    assert hit.gene == "tetA"
    assert hit.database == "tinyamr"
    assert hit.accession == "NC_000913.3:100-900"
    assert hit.function == "TETRACYCLINE"
    assert hit.product == "tetracycline efflux pump TetA"
    assert hit.identity_pct == 100.0
    assert hit.coverage_pct == 100.0
    assert hit.coverage_map == "==============="
    assert hit.strand == "+"


def test_partial_hit_below_default_mincov_is_dropped(tinyamr: Database) -> None:
    """Given a ~50% prefix of blaTEM-1, When screened at mincov=80, Then no hits."""
    report = screen_file(CONTIGS / "partial.fa", tinyamr, PARAMS, dbtype="nucl")
    assert report.hits == ()


def test_partial_hit_kept_at_low_mincov(tinyamr: Database) -> None:
    """Given the same ~50% prefix, When screened at mincov=40, Then the hit
    survives with coverage_pct 50.0 (100*44/88)."""
    params = ScreeningParams(db="tinyamr", mincov=40.0)
    (hit,) = screen_file(CONTIGS / "partial.fa", tinyamr, params, dbtype="nucl").hits
    assert hit.gene == "blaTEM-1"
    assert hit.coverage_pct == 50.0
    assert hit.identity_pct == 100.0


def test_gapped_alignment_reports_broken_map(tinyamr: Database) -> None:
    """Given sul1 with a 3-nt insertion, When screened, Then one hit with
    gap_openings >= 1, gaps=3, coverage 100% ((97-3)/94), and '/' in the map."""
    (hit,) = screen_file(CONTIGS / "gap.fa", tinyamr, PARAMS, dbtype="nucl").hits
    assert hit.gene == "sul1"
    assert hit.gap_openings >= 1
    assert hit.gaps == 3
    assert hit.coverage_pct == 100.0
    assert "/" in hit.coverage_map
    assert len(hit.coverage_map) == 15


def test_unrelated_contig_yields_empty_report(tinyamr: Database) -> None:
    """Given a contig sharing no 11-mer with the db, When screened, Then an
    empty Report with success semantics (hits == ())."""
    report = screen_file(CONTIGS / "none.fa", tinyamr, PARAMS, dbtype="nucl")
    assert report.hits == ()
    assert report.file == str(CONTIGS / "none.fa")


def test_junk_input_raises_input_error(tinyamr: Database, tmp_path: Path) -> None:
    """Given a non-sequence .txt, When screened, Then any2fasta fails and the
    pipeline raises InputError (exit 5)."""
    junk = tmp_path / "junk.txt"
    junk.write_text("this is not sequence data at all\n", encoding="utf-8")
    with pytest.raises(InputError) as excinfo:
        screen_file(junk, tinyamr, PARAMS, dbtype="nucl")
    assert excinfo.value.exit_code == 5
    assert excinfo.value.code == "INVALID_INPUT"


def test_report_is_sorted_by_sequence_then_start(tinyamr: Database) -> None:
    """Given contigB/contigA/contigC with two genes on contigC, When screened,
    Then Report order is (sequence lexicographic, start numeric) — note BLAST
    itself emits contigB, contigA, contigC/sul1, contigC/tetA."""
    report = screen_file(CONTIGS / "sort.fa", tinyamr, PARAMS, dbtype="nucl")
    assert [(hit.sequence, hit.gene) for hit in report.hits] == [
        ("contigA", "tetA"),
        ("contigB", "blaTEM-1"),
        ("contigC", "tetA"),
        ("contigC", "sul1"),
    ]
    assert [hit.start for hit in report.hits] == [1, 1, 1, 90]
