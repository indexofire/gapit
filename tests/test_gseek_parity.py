"""Real-tool parity for the optional gseek engine.

The tests are skipped when gseek or the external binaries are unavailable;
the default pixi environment installs both.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Literal

import pytest

from gapit.blast import run_blastn
from gapit.db import Database, make_blast_db
from gapit.minimap2_run import run_minimap2
from gapit.report import ScreeningParams

pytest.importorskip("gseek")

requires_tools = pytest.mark.skipif(
    any(shutil.which(binary) is None for binary in ("blastn", "makeblastdb", "minimap2")),
    reason="BLAST+ or minimap2 not available",
)

DNA = "ATG" + "ACGT" * 150 + "TAA"
PROTEIN = "M" + "ACDEFGHIKLMNPQRSTVWY" * 25
CODONS = {
    "A": "GCT",
    "C": "TGT",
    "D": "GAT",
    "E": "GAA",
    "F": "TTT",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",
    "L": "CTT",
    "M": "ATG",
    "N": "AAT",
    "P": "CCT",
    "Q": "CAA",
    "R": "CGT",
    "S": "TCT",
    "T": "ACT",
    "V": "GTT",
    "W": "TGG",
    "Y": "TAT",
}


def _db(tmp_path: Path, sequence: str, *, prot: bool = False) -> Database:
    directory = tmp_path / "db"
    directory.mkdir()
    path = directory / "sequences"
    path.write_text(f">gapit|db=test|gene=abc|acc=1|func=test\n{sequence}\n")
    make_blast_db(path, "test", dbtype="prot" if prot else "nucl")
    return Database(name="test", path=directory, sequences_path=path)


def _run(
    query: Path,
    db: Database,
    *,
    dbtype: Literal["nucl", "prot"],
    monkeypatch: pytest.MonkeyPatch,
    mode: str,
):
    monkeypatch.setenv("GAPIT_ENGINE", mode)
    return run_blastn(
        query,
        db,
        ScreeningParams(db="test", minid=80.0, mincov=80.0, threads=1),
        dbtype=dbtype,
    )


@requires_tools
def test_blastn_parity(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    db = _db(tmp_path, DNA)
    query = tmp_path / "q.fa"
    query.write_text(">contig1\n" + "TTTT" + DNA + "GGGG\n")
    legacy = _run(query, db, dbtype="nucl", monkeypatch=monkeypatch, mode="legacy")
    new = _run(query, db, dbtype="nucl", monkeypatch=monkeypatch, mode="gseek")
    assert new == legacy


@requires_tools
def test_blastx_parity(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    db = _db(tmp_path, PROTEIN, prot=True)
    query = tmp_path / "q.fa"
    query.write_text(">contig1\n" + "".join(CODONS[aa] for aa in PROTEIN) + "\n")
    legacy = _run(query, db, dbtype="prot", monkeypatch=monkeypatch, mode="legacy")
    new = _run(query, db, dbtype="prot", monkeypatch=monkeypatch, mode="gseek")
    assert new == legacy


@requires_tools
@pytest.mark.parametrize("read_type", ["sr", "asm20"])
def test_minimap2_parity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, read_type: Literal["sr", "asm20"]
) -> None:
    db = _db(tmp_path, DNA)
    if read_type == "sr":
        query = tmp_path / "r.fq"
        query.write_text("@read1\n" + DNA[30:130] + "\n+\n" + "I" * 100 + "\n")
    else:
        query = tmp_path / "q.fa"
        query.write_text(">contig1\n" + DNA + "\n")
    monkeypatch.setenv("GAPIT_ENGINE", "legacy")
    legacy = run_minimap2([(query, None)], db, read_type=read_type, threads=1, nm_tags=True)
    monkeypatch.setenv("GAPIT_ENGINE", "gseek")
    new = run_minimap2([(query, None)], db, read_type=read_type, threads=1, nm_tags=True)
    # gseek returns a deterministic sorted view; gapit consumers aggregate/sort,
    # so compare the PAF row sets rather than raw emission order.
    assert sorted(new, key=lambda row: row.model_dump_json()) == sorted(
        legacy, key=lambda row: row.model_dump_json()
    )


@requires_tools
def test_minimap2_paired_parity(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    db = _db(tmp_path, DNA)
    r1 = tmp_path / "r1.fq"
    r2 = tmp_path / "r2.fq"
    r1.write_text("@read1\n" + DNA[30:130] + "\n+\n" + "I" * 100 + "\n")
    r2.write_text("@read1\n" + DNA[130:230] + "\n+\n" + "I" * 100 + "\n")
    monkeypatch.setenv("GAPIT_ENGINE", "legacy")
    legacy = run_minimap2([(r1, r2)], db, read_type="sr", threads=1, nm_tags=True)
    monkeypatch.setenv("GAPIT_ENGINE", "gseek")
    new = run_minimap2([(r1, r2)], db, read_type="sr", threads=1, nm_tags=True)
    assert new == legacy
