"""Database layer: discovery, ``~~~`` header parsing, makeblastdb/blastdbcmd wrappers."""

import os
import re
import subprocess
from pathlib import Path
from typing import Literal

from pydantic import BaseModel

from gaita.errors import DatabaseError, DependencyError
from gaita.fasta import iter_fasta

IDSEP = "~~~"


class DbHeader(BaseModel, frozen=True):
    """Parsed ``~~~`` fields of a database sequence id (SPEC.md §4 step 5)."""

    database: str
    gene: str
    accession: str
    resistance: str


class BlastDbInfo(BaseModel, frozen=True):
    """Introspection of one built BLAST database (``blastdbcmd -info``)."""

    n_sequences: int
    dbtype: Literal["nucl", "prot"]
    title: str
    date: str


class Database(BaseModel, frozen=True):
    """One discovered database directory under the datadir."""

    name: str
    path: Path
    sequences_path: Path


class DatabaseInfo(BaseModel, frozen=True):
    """One row of ``gaita list`` output."""

    name: str
    n_sequences: int
    dbtype: Literal["nucl", "prot"]
    date: str


def parse_db_header(seqid: str, default_db: str) -> DbHeader:
    """Split a seqid on ``~~~`` into (database, gene, accession, resistance).

    Without any ``~~~`` the whole seqid is the gene and the default database is
    used; partial headers leave trailing fields empty; an empty database field
    also falls back to ``default_db`` (SPEC.md §4 step 5).
    """
    fields = seqid.split(IDSEP)
    if len(fields) == 1:
        return DbHeader(database=default_db, gene=seqid, accession="", resistance="")
    fields += [""] * (4 - len(fields))
    database, gene, accession, resistance = fields[:4]
    return DbHeader(
        database=database or default_db,
        gene=gene,
        accession=accession,
        resistance=resistance,
    )


_AGTC_TABLE = str.maketrans("", "", "AGTCagtc")


def mol_type(letters: str) -> Literal["nucl", "prot"]:
    """Abricate heuristic (SPEC.md §2): protein iff non-[AGTC] letters are
    strictly more than half of the input (empty input is nucleotide)."""
    non_agtc = len(letters.translate(_AGTC_TABLE))
    return "prot" if 2 * non_agtc > len(letters) else "nucl"


def _run(argv: list[str]) -> subprocess.CompletedProcess[str]:
    """Run an external tool with an argv list (never a shell)."""
    try:
        return subprocess.run(argv, check=False, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise DependencyError(
            f"required binary not found on PATH: {argv[0]}",
            code="MISSING_DEPENDENCY",
        ) from exc


def make_blast_db(sequences_path: Path, name: str) -> None:
    """(Re)build the BLAST index for one database directory."""
    letters = "".join(record.sequence for record in iter_fasta(sequences_path))
    dbtype = mol_type(letters)
    for index_file in sequences_path.parent.glob(f"{sequences_path.name}.[np]??"):
        index_file.unlink()
    result = _run(
        [
            "makeblastdb",
            "-in",
            str(sequences_path),
            "-title",
            name,
            "-dbtype",
            dbtype,
            "-logfile",
            "/dev/null",
        ]
    )
    if result.returncode != 0:
        raise DatabaseError(
            f"makeblastdb failed for {name}: {result.stderr.strip()}",
            code="MAKEBLASTDB_FAILED",
        )


_SEQUENCES_RE = re.compile(r"([\d,]+)\s+sequences;")
_TITLE_RE = re.compile(r"^Database:[ \t]*(.*)$", re.MULTILINE)
_DATE_RE = re.compile(r"Date:\s*([A-Za-z]{3})\s+(\d{1,2}),\s+(\d{4})")


def parse_blastdbcmd_info(text: str) -> BlastDbInfo:
    """Parse ``blastdbcmd -info`` output (SPEC.md §2)."""
    sequences_match = _SEQUENCES_RE.search(text)
    title_match = _TITLE_RE.search(text)
    date_match = _DATE_RE.search(text)
    if sequences_match is None or title_match is None or date_match is None:
        raise DatabaseError(
            f"could not parse blastdbcmd output: {text!r}",
            code="BLASTDBCMD_PARSE_FAILED",
        )
    month, day, year = date_match.group(1), date_match.group(2), date_match.group(3)
    return BlastDbInfo(
        n_sequences=int(sequences_match.group(1).replace(",", "")),
        dbtype="prot" if "total residues" in text else "nucl",
        title=title_match.group(1).strip(),
        date=f"{year}-{month}-{int(day):02d}",
    )


def blast_db_info(db_prefix: Path) -> BlastDbInfo:
    """Introspect a built BLAST database via ``blastdbcmd -info``."""
    result = _run(["blastdbcmd", "-info", "-db", str(db_prefix)])
    if result.returncode != 0:
        raise DatabaseError(
            f"Database {db_prefix} is not indexed, please try: gaita setupdb",
            code="DATABASE_NOT_INDEXED",
        )
    return parse_blastdbcmd_info(result.stdout)


def discover_databases(datadir: Path) -> list[Database]:
    """Immediate subdirectories of the datadir holding a readable ``sequences``
    file, sorted by name."""
    databases: list[Database] = []
    for child in datadir.iterdir():
        sequences = child / "sequences"
        if child.is_dir() and sequences.is_file() and os.access(sequences, os.R_OK):
            databases.append(Database(name=child.name, path=child, sequences_path=sequences))
    return sorted(databases, key=lambda database: database.name)


def list_databases(datadir: Path, *, setupdb: bool) -> list[DatabaseInfo]:
    """Enumerate databases (building indices first when ``setupdb``), requiring
    every database to be indexed; returns one DatabaseInfo per database."""
    infos: list[DatabaseInfo] = []
    for database in discover_databases(datadir):
        if setupdb:
            make_blast_db(database.sequences_path, database.name)
        sequences = database.sequences_path
        index_exists = any(
            (sequences.parent / f"{sequences.name}{suffix}").exists() for suffix in (".nin", ".pin")
        )
        if not index_exists:
            raise DatabaseError(
                f"Database {database.name} is not indexed, please try: gaita setupdb",
                code="DATABASE_NOT_INDEXED",
            )
        info = blast_db_info(sequences)
        infos.append(
            DatabaseInfo(
                name=database.name,
                n_sequences=info.n_sequences,
                dbtype=info.dbtype,
                date=info.date,
            )
        )
    return infos
