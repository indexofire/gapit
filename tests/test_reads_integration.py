"""Integration tests: real minimap2 read screening on committed fixtures.

The tinyreads db carries longer genes (tetX 522 nt, sulY 261 nt) built from
tinyamr material because minimap2 -x sr systematically soft-clips ~5 nt at
each end of alignments on the 79-94 nt tinyamr genes, capping breadth below
the 90% presence threshold (validated empirically; see phase notes).
"""

import gzip
import shutil
from pathlib import Path

import pytest

from gapit.db import Database, discover_databases
from gapit.errors import InputError
from gapit.reads import ReadFileKind, detect_read_kind, screen_reads

READS_DB_DIR = Path(__file__).parent / "data" / "reads_db"
TINYAMR_DB_DIR = Path(__file__).parent / "data" / "db"
READS = Path(__file__).parent / "data" / "reads"


@pytest.fixture()
def readsdb(tmp_path: Path) -> Database:
    datadir = tmp_path / "datadir"
    shutil.copytree(READS_DB_DIR, datadir)
    (database,) = discover_databases(datadir)
    return database


@pytest.fixture()
def tinyamr(tmp_path: Path) -> Database:
    datadir = tmp_path / "datadir"
    shutil.copytree(TINYAMR_DB_DIR, datadir)
    (database,) = discover_databases(datadir)
    return database


def test_full_coverage_gene_is_present(readsdb: Database) -> None:
    """Given tiled 100 nt reads across tetX, When screened with -x sr, Then
    tetX is present with near-total breadth and all 12 reads mapped."""
    report = screen_reads(
        [(READS / "tetx_full.fq", None)], readsdb, read_type="sr", min_breadth=90.0, threads=1
    )
    (entry,) = report.genes
    assert entry.gene == "tetX"
    assert entry.database == "tinyreads"
    assert entry.accession == "SYN-001"
    assert entry.function == "TETRACYCLINE"
    assert entry.product == "extended resistance determinant tetX"
    assert entry.present is True
    assert entry.breadth_pct >= 95.0
    assert entry.reads_mapped == 12
    assert entry.tlen == 522
    assert report.reads == (str(READS / "tetx_full.fq"),)


def test_partial_coverage_gene_is_absent(readsdb: Database) -> None:
    """Given reads covering only ~60% of sulY, When screened, Then sulY is
    listed but present is False."""
    report = screen_reads(
        [(READS / "suly_partial.fq", None)], readsdb, read_type="sr", min_breadth=90.0, threads=1
    )
    (entry,) = report.genes
    assert entry.gene == "sulY"
    assert entry.present is False
    assert 0.0 < entry.breadth_pct < 90.0


def test_junk_reads_detect_no_genes(readsdb: Database) -> None:
    """Given unrelated reads, When screened, Then zero genes and success."""
    report = screen_reads(
        [(READS / "junk.fq", None)], readsdb, read_type="sr", min_breadth=90.0, threads=1
    )
    assert report.genes == ()


def test_gzipped_reads(readsdb: Database) -> None:
    """Given the same reads gzip-compressed, When screened, Then identical
    presence behavior (minimap2 reads .gz natively)."""
    report = screen_reads(
        [(READS / "tetx_full.fq.gz", None)], readsdb, read_type="sr", min_breadth=90.0, threads=1
    )
    (entry,) = report.genes
    assert entry.gene == "tetX"
    assert entry.present is True
    assert entry.reads_mapped == 12


def test_paired_lane(readsdb: Database) -> None:
    """Given the tiled reads split into two mate files, When screened as one
    paired lane, Then the union coverage and read count match the single-file
    run, and report.reads lists r1 then r2."""
    report = screen_reads(
        [(READS / "tetx_R1.fq", READS / "tetx_R2.fq")],
        readsdb,
        read_type="sr",
        min_breadth=90.0,
        threads=1,
    )
    (entry,) = report.genes
    assert entry.present is True
    assert entry.reads_mapped == 12
    assert report.reads == (str(READS / "tetx_R1.fq"), str(READS / "tetx_R2.fq"))


def test_lane_union_crosses_threshold(readsdb: Database) -> None:
    """Given two lanes each covering only ~half of tetX, When screened
    separately, Then each is absent (~50% breadth); When screened together as
    one sample, Then the union (~97.7%) crosses the 90% threshold."""
    for lane in ("tetx_lane1.fq", "tetx_lane2.fq"):
        solo = screen_reads(
            [(READS / lane, None)], readsdb, read_type="sr", min_breadth=90.0, threads=1
        )
        (entry,) = solo.genes
        assert entry.present is False
        assert entry.breadth_pct < 90.0
    union = screen_reads(
        [(READS / "tetx_lane1.fq", None), (READS / "tetx_lane2.fq", None)],
        readsdb,
        read_type="sr",
        min_breadth=90.0,
        threads=1,
    )
    (entry,) = union.genes
    assert entry.gene == "tetX"
    assert entry.present is True
    assert entry.breadth_pct >= 95.0
    assert entry.reads_mapped == 12
    assert union.reads == (str(READS / "tetx_lane1.fq"), str(READS / "tetx_lane2.fq"))


def test_paired_lane_union(readsdb: Database) -> None:
    """Given two paired lanes (2x2 mate files), When screened, Then the gene is
    found once with reads from all four files counted."""
    report = screen_reads(
        [
            (READS / "tetx_pe1_R1.fq", READS / "tetx_pe1_R2.fq"),
            (READS / "tetx_pe2_R1.fq", READS / "tetx_pe2_R2.fq"),
        ],
        readsdb,
        read_type="sr",
        min_breadth=90.0,
        threads=1,
    )
    (entry,) = report.genes
    assert entry.present is True
    assert entry.breadth_pct >= 95.0
    assert entry.reads_mapped == 12
    assert len(report.reads) == 4


def test_tinyamr_partial_below_threshold(tinyamr: Database) -> None:
    """Given a 44 nt read on blaTEM-1 (88 nt), When screened against tinyamr,
    Then blaTEM-1 is listed with present False."""
    report = screen_reads(
        [(READS / "bla_partial.fq", None)], tinyamr, read_type="sr", min_breadth=90.0, threads=1
    )
    (entry,) = report.genes
    assert entry.gene == "blaTEM-1"
    assert entry.present is False
    assert entry.breadth_pct < 90.0


def test_detect_read_kind_fastq(tmp_path: Path) -> None:
    """Given a file whose first byte is '@', When detected, Then FASTQ."""
    path = tmp_path / "r1.fq"
    path.write_text("@read1\nACGT\n+\nIIII\n", encoding="utf-8")
    assert detect_read_kind(path) is ReadFileKind.fastq


def test_detect_read_kind_fasta(tmp_path: Path) -> None:
    """Given a file whose first byte is '>', When detected, Then FASTA."""
    path = tmp_path / "assembly.fa"
    path.write_text(">contig1\nACGT\n", encoding="utf-8")
    assert detect_read_kind(path) is ReadFileKind.fasta


def test_detect_read_kind_skips_leading_whitespace(tmp_path: Path) -> None:
    """Given a FASTQ preceded by blank lines and spaces, When detected, Then
    FASTQ (peek walks to the first non-whitespace byte)."""
    path = tmp_path / "padded.fq"
    path.write_text("\n\n   @read1\nACGT\n+\nIIII\n", encoding="utf-8")
    assert detect_read_kind(path) is ReadFileKind.fastq


def test_detect_read_kind_gzipped_fastq(tmp_path: Path) -> None:
    """Given a gzip-wrapped FASTQ (magic 1f 8b first), When detected, Then
    FASTQ — the peek goes through the decompressor, matching minimap2's
    native .gz support."""
    path = tmp_path / "r1.fq.gz"
    with gzip.open(path, "wb") as handle:
        handle.write(b"@read1\nACGT\n+\nIIII\n")
    assert detect_read_kind(path) is ReadFileKind.fastq


def test_detect_read_kind_gzipped_fasta(tmp_path: Path) -> None:
    """Given a gzip-wrapped FASTA, When detected, Then FASTA."""
    path = tmp_path / "assembly.fa.gz"
    with gzip.open(path, "wb") as handle:
        handle.write(b">contig1\nACGT\n")
    assert detect_read_kind(path) is ReadFileKind.fasta


def test_detect_read_kind_rejects_garbage(tmp_path: Path) -> None:
    """Given a file whose first non-whitespace byte is neither '>' nor '@',
    When detected, Then typed InputError with the file in context."""
    path = tmp_path / "garbage.txt"
    path.write_text("Nonsense, not sequencing data\n", encoding="utf-8")
    with pytest.raises(InputError) as excinfo:
        detect_read_kind(path)
    assert excinfo.value.code == "INVALID_READS_FORMAT"
    assert excinfo.value.context["file"] == str(path)


def test_detect_read_kind_rejects_empty(tmp_path: Path) -> None:
    """Given an empty (or whitespace-only) file, When detected, Then typed
    InputError."""
    path = tmp_path / "empty.fq"
    path.write_text("   \n", encoding="utf-8")
    with pytest.raises(InputError) as excinfo:
        detect_read_kind(path)
    assert excinfo.value.code == "INVALID_READS_FORMAT"
