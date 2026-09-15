"""Unit tests for database discovery and blastdbcmd output parsing (gaita.db)."""

from pathlib import Path

import pytest

from gaita.db import discover_databases, list_databases, parse_blastdbcmd_info
from gaita.errors import DatabaseError

# Captured verbatim from `blastdbcmd -info -db <fixture>/sequences` (BLAST+ 2.17,
# tab-indented, "Longest sequence" sharing the Date line is authentic).
NUCL_SAMPLE = (
    "Database: tinyamr\n"
    "\t3 sequences; 261 total bases\n"
    "\n"
    "Date: Sep 15, 2026  11:42 PM\tLongest sequence: 94 bases\n"
    "\n"
    "BLASTDB Version: 5\n"
    "\n"
    "Volumes:\n"
    "\t/tmp/opencode/nucl/sequences\n"
)

# Same shape as a real protein db (`total residues` instead of `total bases`),
# with thousands separators and the double-space day variant seen from BLAST+.
PROT_SAMPLE = (
    "Database: Protein AMR db\n"
    "\t12,345 sequences; 1,234,567 total residues\n"
    "\n"
    "Date: Jan  5, 2026  10:00:00\tLongest sequence: 3,654 residues\n"
    "\n"
    "BLASTDB Version: 5\n"
    "\n"
    "Volumes:\n"
    "\t/tmp/opencode/prot/sequences\n"
)


def test_parse_nucl_sample() -> None:
    info = parse_blastdbcmd_info(NUCL_SAMPLE)
    assert info.n_sequences == 3
    assert info.dbtype == "nucl"
    assert info.title == "tinyamr"
    assert info.date == "2026-Sep-15"


def test_parse_prot_sample() -> None:
    info = parse_blastdbcmd_info(PROT_SAMPLE)
    assert info.n_sequences == 12345
    assert info.dbtype == "prot"
    assert info.title == "Protein AMR db"
    assert info.date == "2026-Jan-05"


def write_sequences(db_dir: Path, content: str = ">s\nACGT\n") -> Path:
    db_dir.mkdir(parents=True)
    (db_dir / "sequences").write_text(content, encoding="utf-8")
    return db_dir


def test_discover_databases_finds_only_dirs_with_sequences_file(tmp_path: Path) -> None:
    """Given a datadir mixing real dbs, an empty dir, a dir whose `sequences` is a
    directory, and a plain file, When discovered, Then only real dbs remain."""
    aaa = write_sequences(tmp_path / "aaa")
    (tmp_path / "bbb").mkdir()  # no sequences file
    (tmp_path / "ccc").mkdir()
    (tmp_path / "ccc" / "sequences").mkdir()  # sequences is itself a directory
    (tmp_path / "plainfile").write_text("not a dir", encoding="utf-8")
    dbs = discover_databases(tmp_path)
    assert [db.name for db in dbs] == ["aaa"]
    assert dbs[0].path == aaa
    assert dbs[0].sequences_path == aaa / "sequences"


def test_discover_databases_sorted_by_name(tmp_path: Path) -> None:
    for name in ("zeta", "alpha", "mid"):
        write_sequences(tmp_path / name)
    assert [db.name for db in discover_databases(tmp_path)] == ["alpha", "mid", "zeta"]


def test_list_databases_without_setupdb_raises_not_indexed(tmp_path: Path) -> None:
    """Given a discovered db without a BLAST index, When listed, Then DatabaseError
    (exit 4) pointing at `gaita setupdb`."""
    write_sequences(tmp_path / "tinyamr")
    with pytest.raises(DatabaseError) as excinfo:
        list_databases(tmp_path, setupdb=False)
    assert "not indexed" in str(excinfo.value)
    assert "gaita setupdb" in str(excinfo.value)
    assert excinfo.value.exit_code == 4
