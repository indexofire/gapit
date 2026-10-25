"""Unit tests for the optional gseek bridge and engine dispatch."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import pytest

from gapit import blast, minimap2_run
from gapit.db import Database
from gapit.gseek_engine import (
    engine_mode,
    iter_paf_gseek,
    run_blast_gseek,
)
from gapit.paf import PafRecord
from gapit.report import ScreeningParams


@dataclass
class _Record:
    query_id: str = "contig1"
    target_id: str = "gapit|db=test|gene=abc|acc=1|func=x"
    query_length: int | None = 100
    target_length: int | None = 100
    query_start: int | None = 10
    query_end: int | None = 20
    target_start: int | None = 30
    target_end: int | None = 40
    identity: float | None = 99.0
    aligned_length: int | None = 11
    matches: int | None = 10
    gap_count: int | None = 1
    gap_openings: int | None = 1
    evalue: float | None = 1e-20
    mapq: int | None = 42
    target_description: str | None = "abc product"
    is_primary: bool = True
    tags: dict[str, str] | None = None
    reverse: bool = False

    def __post_init__(self) -> None:
        if self.tags is None:
            self.tags = {}

    @property
    def is_reverse(self) -> bool:
        return self.reverse


class _FakeGseek:
    def __init__(self, records: list[_Record]) -> None:
        self.records = records
        self.kwargs: dict[str, object] = {}

    def align(self, **kwargs: object) -> SimpleNamespace:
        self.kwargs = kwargs
        return SimpleNamespace(records=tuple(self.records))


def _database(tmp_path: Path) -> Database:
    return Database(name="test", path=tmp_path, sequences_path=tmp_path / "sequences")


def test_engine_mode_default_and_validation(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("GAPIT_ENGINE", raising=False)
    assert engine_mode() == "legacy"
    monkeypatch.setenv("GAPIT_ENGINE", "gseek")
    assert engine_mode() == "gseek"
    monkeypatch.setenv("GAPIT_ENGINE", "bogus")
    with pytest.raises(Exception, match="unknown engine mode"):
        engine_mode()


def test_blast_mapping_preserves_blast_coordinates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = _FakeGseek([_Record()])
    monkeypatch.setattr("gapit.gseek_engine._gseek", lambda: fake)
    rows = run_blast_gseek(
        ">contig1\nACGT\n",
        _database(tmp_path),
        ScreeningParams(db="test", threads=2),
        dbtype="nucl",
    )
    assert len(rows) == 1
    row = rows[0]
    assert (row.qseqid, row.qstart, row.qend, row.qlen) == ("contig1", 11, 20, 100)
    assert (row.sseqid, row.sstart, row.send, row.slen) == (
        "gapit|db=test|gene=abc|acc=1|func=x",
        31,
        40,
        100,
    )
    assert row.sstrand == "plus"
    assert (row.evalue, row.length, row.pident, row.gaps, row.gapopen) == (
        1e-20,
        11,
        99.0,
        1,
        1,
    )
    assert row.stitle == "abc product"
    assert fake.kwargs["threads"] == 2
    assert fake.kwargs["record_manifest"] is False


def test_blast_mapping_reverse_and_protein_params(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fake = _FakeGseek([_Record(reverse=True)])
    monkeypatch.setattr("gapit.gseek_engine._gseek", lambda: fake)
    rows = run_blast_gseek(
        ">contig1\nACGT\n",
        _database(tmp_path),
        ScreeningParams(db="test"),
        dbtype="prot",
    )
    assert (rows[0].sstart, rows[0].send, rows[0].sstrand) == (40, 31, "N/A")
    params = fake.kwargs["params"]
    assert isinstance(params, dict)
    assert params["task"] == "blastx-fast"
    assert params["seg"] is False
    assert "min_identity" not in params


def test_paf_mapping_and_paired_query(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    record = _Record(
        query_id="read1",
        target_id="gene1",
        query_start=0,
        query_end=99,
        target_start=5,
        target_end=104,
        matches=99,
        aligned_length=100,
        mapq=7,
        tags={"NM": "1", "cs": ":100"},
    )
    fake = _FakeGseek([record])
    monkeypatch.setattr("gapit.gseek_engine._gseek", lambda: fake)
    r1 = tmp_path / "r1.fq"
    r2 = tmp_path / "r2.fq"
    r1.write_text("@r1\nACGT\n+\nIIII\n")
    r2.write_text("@r2\nACGT\n+\nIIII\n")
    rows = list(
        iter_paf_gseek(
            [(r1, r2)],
            _database(tmp_path),
            read_type="sr",
            threads=2,
            nm_tags=True,
        )
    )
    assert rows == [
        PafRecord(
            qname="read1",
            qlen=100,
            qstart=0,
            qend=99,
            strand="+",
            tname="gene1",
            tlen=100,
            tstart=5,
            tend=104,
            nmatch=99,
            alen=100,
            mapq=7,
            is_primary=True,
            nm=1,
            cs=":100",
        )
    ]
    assert fake.kwargs["query"] == (r1, r2)
    assert fake.kwargs["paired"] is True
    params = fake.kwargs["params"]
    assert isinstance(params, dict)
    assert params == {"preset": "sr", "secondary": "all", "cs": True, "k": 11, "w": 6, "B": 2}


def _blast_line() -> str:
    return "\t".join(
        [
            "contig1",
            "1",
            "10",
            "100",
            "gapit|db=test|gene=abc|acc=1|func=x",
            "1",
            "10",
            "100",
            "plus",
            "1e-20",
            "10",
            "100.0",
            "0",
            "0",
            "abc product",
        ]
    )


def test_run_blastn_dispatches_to_gseek(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    query = tmp_path / "q.fa"
    query.write_text(">contig1\nACGT\n")
    database = _database(tmp_path)
    params = ScreeningParams(db="test")
    sentinel = blast.BlastRow(
        qseqid="contig1",
        qstart=1,
        qend=10,
        qlen=100,
        sseqid="gene1",
        sstart=1,
        send=10,
        slen=100,
        sstrand="plus",
        evalue=1e-20,
        length=10,
        pident=100.0,
        gaps=0,
        gapopen=0,
        stitle="x",
    )
    monkeypatch.setenv("GAPIT_ENGINE", "gseek")

    def _new(*args: object, **kwargs: object) -> list[blast.BlastRow]:
        return [sentinel]

    def _legacy(*args: object, **kwargs: object) -> str:
        pytest.fail("legacy path ran")

    monkeypatch.setattr(blast, "run_blast_gseek", _new)
    monkeypatch.setattr(blast, "_pipeline", _legacy)
    assert blast.run_blastn(query, database, params, dbtype="nucl") == [sentinel]


def test_run_blastn_shadow_returns_legacy_and_logs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    query = tmp_path / "q.fa"
    query.write_text(">contig1\nACGT\n")
    database = _database(tmp_path)
    params = ScreeningParams(db="test")
    monkeypatch.setenv("GAPIT_ENGINE", "shadow")

    def _legacy(*args: object, **kwargs: object) -> str:
        return _blast_line()

    monkeypatch.setattr(blast, "_pipeline", _legacy)
    new = blast.BlastRow(
        qseqid="contig1",
        qstart=1,
        qend=10,
        qlen=100,
        sseqid="OTHER",
        sstart=1,
        send=10,
        slen=100,
        sstrand="plus",
        evalue=1e-20,
        length=10,
        pident=90.0,
        gaps=0,
        gapopen=0,
        stitle="other",
    )

    def _new(*args: object, **kwargs: object) -> list[blast.BlastRow]:
        return [new]

    monkeypatch.setattr(blast, "run_blast_gseek", _new)
    with caplog.at_level(logging.WARNING, logger="gapit.gseek"):
        result = blast.run_blastn(query, database, params, dbtype="nucl")
    assert result[0].sseqid.startswith("gapit|")
    assert "legacy and gseek results differ" in caplog.text


def test_iter_minimap2_dispatches_to_gseek(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    database = _database(tmp_path)
    expected = PafRecord(
        qname="read1",
        qlen=100,
        qstart=0,
        qend=99,
        strand="+",
        tname="gene1",
        tlen=100,
        tstart=1,
        tend=100,
        nmatch=99,
        alen=100,
        mapq=60,
    )
    monkeypatch.setenv("GAPIT_ENGINE", "gseek")

    def _rows(*args: object, **kwargs: object) -> object:
        return iter([expected])

    monkeypatch.setattr(minimap2_run, "iter_paf_gseek", _rows)
    rows = list(
        minimap2_run.iter_minimap2([(tmp_path / "r.fq", None)], database, read_type="sr", threads=1)
    )
    assert rows == [expected]


def test_iter_minimap2_shadow_returns_legacy_and_logs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    database = _database(tmp_path)
    legacy = PafRecord(
        qname="read1",
        qlen=100,
        qstart=0,
        qend=99,
        strand="+",
        tname="gene1",
        tlen=100,
        tstart=1,
        tend=100,
        nmatch=99,
        alen=100,
        mapq=60,
    )
    new = legacy.model_copy(update={"mapq": 0})
    monkeypatch.setenv("GAPIT_ENGINE", "shadow")

    def _legacy(*args: object, **kwargs: object) -> object:
        return iter([legacy])

    monkeypatch.setattr(minimap2_run, "_iter_legacy", _legacy)

    def _rows(*args: object, **kwargs: object) -> object:
        return iter([new])

    monkeypatch.setattr(minimap2_run, "iter_paf_gseek", _rows)
    with caplog.at_level(logging.WARNING, logger="gapit.gseek"):
        rows = list(
            minimap2_run.iter_minimap2(
                [(tmp_path / "r.fq", None)], database, read_type="sr", threads=1
            )
        )
    assert rows == [legacy]
    assert "legacy and gseek results differ" in caplog.text


def test_run_blastn_shadow_falls_back_on_gseek_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    query = tmp_path / "q.fa"
    query.write_text(">contig1\nACGT\n")
    database = _database(tmp_path)
    params = ScreeningParams(db="test")
    monkeypatch.setenv("GAPIT_ENGINE", "shadow")

    def _legacy(*args: object, **kwargs: object) -> str:
        return _blast_line()

    def _fail(*args: object, **kwargs: object) -> object:
        raise blast.GapitError("boom", code="GSEEK_FAILED")

    monkeypatch.setattr(blast, "_pipeline", _legacy)
    monkeypatch.setattr(blast, "run_blast_gseek", _fail)
    with caplog.at_level(logging.WARNING, logger="gapit.gseek"):
        rows = blast.run_blastn(query, database, params, dbtype="nucl")
    assert rows[0].sseqid.startswith("gapit|")
    assert "gseek path failed" in caplog.text


def test_iter_minimap2_shadow_falls_back_on_gseek_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    database = _database(tmp_path)
    expected = PafRecord(
        qname="read1",
        qlen=100,
        qstart=0,
        qend=99,
        strand="+",
        tname="gene1",
        tlen=100,
        tstart=1,
        tend=100,
        nmatch=99,
        alen=100,
        mapq=60,
    )
    monkeypatch.setenv("GAPIT_ENGINE", "shadow")

    def _legacy(*args: object, **kwargs: object) -> object:
        return iter([expected])

    def _fail(*args: object, **kwargs: object) -> object:
        raise minimap2_run.GapitError("boom", code="GSEEK_FAILED")

    monkeypatch.setattr(minimap2_run, "_iter_legacy", _legacy)
    monkeypatch.setattr(minimap2_run, "iter_paf_gseek", _fail)
    with caplog.at_level(logging.WARNING, logger="gapit.gseek"):
        rows = list(
            minimap2_run.iter_minimap2(
                [(tmp_path / "r.fq", None)], database, read_type="sr", threads=1
            )
        )
    assert rows == [expected]
    assert "gseek path failed" in caplog.text


def test_minimap2_long_target_names_disable_index(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    long = Database(name="long", path=tmp_path, sequences_path=tmp_path / "long.fa")
    long.sequences_path.write_text(">" + "x" * 256 + "\nACGT\n")
    fake = _FakeGseek([])
    monkeypatch.setattr("gapit.gseek_engine._gseek", lambda: fake)
    rows = list(
        iter_paf_gseek(
            [(tmp_path / "r.fq", None)],
            long,
            read_type="sr",
            threads=1,
        )
    )
    assert rows == []
    assert fake.kwargs["index"] == "none"
