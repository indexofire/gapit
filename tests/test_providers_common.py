"""Tests for the Wave B0 shared provider infrastructure (providers/common.py).

fetch_provider runs the generic abricate-get_db pipeline end to end offline:
a synthetic Provider with a file:// source URL exercises download -> transform
-> normalize -> dedupe -> sort -> records.jsonl -> build_database against the
pixi env's real makeblastdb/minimap2 (same pattern as test_dbbuild.py).

Private helpers (_download, _normalize_*, _dedupe) are covered through the
public API: basedpyright strict (reportPrivateUsage) blocks tests from
importing underscore members, the Wave A3 lesson."""

import urllib.request
from collections.abc import Iterable
from pathlib import Path
from typing import Literal, Self

import pytest

from gapit.errors import DatabaseError
from gapit.fasta import iter_fasta
from gapit.providers.common import Provider, fetch_provider
from gapit.records import Record, read_manifest, read_records

SYN = "synamr"
FETCHED_AT = "2026-09-17T12:00:00Z"

# 60 bp sequences (Wave A3 pinned minimap2 indexing 60-70 bp synthetic seqs).
SEQ_DUP = "ACGT" * 15
SEQ_NORM_RAW = "acgtACGT-nnN" * 5  # lowercase + non-AGCT junk on purpose
SEQ_NORM = "ACGTACGTNNNN" * 5
SEQ_PLAIN = "ACGTAG" * 10

UPSTREAM_FASTA = (
    f">zzz_dup first occurrence wins\n{SEQ_DUP}\n"
    f">aaa_dup duplicate sequence is dropped\n{SEQ_DUP}\n"
    f">mmm_norm junk letters get normalized\n{SEQ_NORM_RAW}\n"
    f">bbb_plain boring sequence\n{SEQ_PLAIN}\n"
)

PROT_FASTA = f">ppp_x\n{'MkvGly-28' * 6}\n>aaa_y\n{'MWVqea' * 8}\n"
SEQ_PROT_NORM = "MKVGLYXXX" * 6  # '-' '2' '8' -> X after uppercase

# Transform-order metadata: function classes arrive unsorted with whitespace runs.
FUNCTION: dict[str, tuple[str, ...]] = {
    "zzz_dup": ("zeta  class", "alpha"),  # sorted -> alpha first; '  ' -> '_'
    "bbb_plain": ("single",),
}


def syn_transform(workdir: Path) -> Iterable[Record]:
    """Load the downloaded upstream fasta and yield Records (db = SYN)."""
    for fasta in iter_fasta(workdir / "upstream.fa"):
        yield Record(
            db=SYN,
            gene=fasta.id,
            sequence=fasta.sequence,
            function=FUNCTION.get(fasta.id, ()),
        )


def empty_transform(workdir: Path) -> Iterable[Record]:
    """A provider transform that legitimately yields nothing."""
    return ()


def write_source(tmp_path: Path, content: str = UPSTREAM_FASTA) -> str:
    """Write an upstream fasta fixture into tmp_path and return its file:// URL."""
    source = tmp_path / "upstream.fa"
    source.write_text(content, encoding="utf-8")
    return source.as_uri()


def syn_provider(url: str, dbtype: Literal["nucl", "prot"] = "nucl") -> Provider:
    return Provider(
        name=SYN,
        description="synthetic Wave B0 test provider",
        source_urls=(url,),
        dbtype=dbtype,
        transform=syn_transform,
    )


def test_fetch_provider_runs_full_pipeline(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Given a synthetic file:// provider whose transform yields 4 records
    (one duplicate sequence, one junk-laden sequence, unsorted function
    classes, out-of-order genes), When fetched, Then every artifact is built,
    records are deduped-first-wins / normalized / sorted by gene, the manifest
    matches, the workdir is cleaned, and progress went to stderr only."""
    url = write_source(tmp_path)
    db_dir = tmp_path / "datadir" / SYN  # does not exist yet: fetch creates it
    manifest = fetch_provider(syn_provider(url), db_dir, fetched_at=FETCHED_AT, quiet=False)

    for name in (
        "records.jsonl",
        "sequences",
        "sequences.nin",
        "sequences.mmi",
        "gapit-manifest.json",
    ):
        assert (db_dir / name).is_file(), name

    assert manifest.name == SYN
    assert manifest.n_records == 3
    assert manifest.dbtype == "nucl"
    assert manifest.source_urls == (url,)
    assert manifest.fetched_at == FETCHED_AT
    assert read_manifest(db_dir / "gapit-manifest.json") == manifest

    records = tuple(read_records(db_dir / "records.jsonl"))
    assert [record.gene for record in records] == [
        "bbb_plain",
        "mmm_norm",
        "zzz_dup",
    ]
    by_gene = {record.gene: record for record in records}
    # duplicate SEQUENCE dropped and the FIRST (lexically larger!) one won
    assert "aaa_dup" not in by_gene
    assert by_gene["zzz_dup"].sequence == SEQ_DUP
    # function classes sorted, whitespace runs collapsed to single '_'
    assert by_gene["zzz_dup"].function == ("alpha", "zeta_class")
    # sequence uppercased, every non-AGCT letter -> N
    assert by_gene["mmm_norm"].sequence == SEQ_NORM
    assert all(record.db == SYN for record in records)

    # workdir cleaned: every entry left in db_dir is a built artifact file
    assert all(entry.is_file() for entry in db_dir.iterdir())

    captured = capsys.readouterr()
    assert captured.out == ""  # stdout purity
    assert "dropped 1 duplicate" in captured.err
    assert "kept 3" in captured.err


def test_fetch_provider_refuses_overwrite_without_force(tmp_path: Path) -> None:
    """Given an already-built database, When fetched again without force,
    Then DB_ALREADY_EXISTS (exit 4) names the db — parity with upstream
    'Won't overwrite existing (use --force)'."""
    url = write_source(tmp_path)
    db_dir = tmp_path / SYN
    db_dir.mkdir()
    fetch_provider(syn_provider(url), db_dir, fetched_at=FETCHED_AT)
    with pytest.raises(DatabaseError) as excinfo:
        fetch_provider(syn_provider(url), db_dir, fetched_at=FETCHED_AT)
    error = excinfo.value
    assert error.code == "DB_ALREADY_EXISTS"
    assert error.exit_code == 4
    assert error.context == {"db": SYN}


def test_fetch_provider_force_rebuilds_existing_database(tmp_path: Path) -> None:
    """Given an already-built database, When fetched with force=True, Then the
    rebuild succeeds and the artifact set is complete again."""
    url = write_source(tmp_path)
    db_dir = tmp_path / SYN
    db_dir.mkdir()
    fetch_provider(syn_provider(url), db_dir, fetched_at=FETCHED_AT)
    manifest = fetch_provider(syn_provider(url), db_dir, fetched_at=FETCHED_AT, force=True)
    assert manifest.n_records == 3
    assert read_manifest(db_dir / "gapit-manifest.json") == manifest
    assert (db_dir / "sequences.mmi").is_file()


def test_fetch_provider_empty_transform_raises(tmp_path: Path) -> None:
    """Given a provider whose transform yields zero records, When fetched,
    Then PROVIDER_EMPTY (exit 4) fires and neither records.jsonl nor a
    manifest was written."""
    url = write_source(tmp_path)
    provider = Provider(
        name=SYN,
        description="synthetic empty provider",
        source_urls=(url,),
        dbtype="nucl",
        transform=empty_transform,
    )
    db_dir = tmp_path / SYN
    with pytest.raises(DatabaseError) as excinfo:
        fetch_provider(provider, db_dir, fetched_at=FETCHED_AT)
    error = excinfo.value
    assert error.code == "PROVIDER_EMPTY"
    assert error.exit_code == 4
    assert error.context == {"db": SYN}
    assert not (db_dir / "records.jsonl").exists()
    assert not (db_dir / "gapit-manifest.json").exists()


def test_fetch_provider_missing_source_raises_download_failed(
    tmp_path: Path,
) -> None:
    """Given a source URL that does not exist, When fetched, Then
    DOWNLOAD_FAILED (exit 4) carries the URL in context and the workdir is
    cleaned up (no directories left in db_dir)."""
    url = (tmp_path / "no-such-file.fa").as_uri()
    db_dir = tmp_path / SYN
    with pytest.raises(DatabaseError) as excinfo:
        fetch_provider(syn_provider(url), db_dir, fetched_at=FETCHED_AT)
    error = excinfo.value
    assert error.code == "DOWNLOAD_FAILED"
    assert error.exit_code == 4
    assert error.context == {"url": url}
    assert not (db_dir / "gapit-manifest.json").exists()
    assert all(entry.is_file() for entry in db_dir.iterdir())


def test_fetch_provider_prot_normalizes_to_x_and_skips_mmi(
    tmp_path: Path,
) -> None:
    """Given a prot provider whose sequences carry lowercase and non-A-Z
    characters, When fetched, Then sequences are uppercased with non-A-Z ->
    X, a .pin index exists, no .mmi is created, and the manifest records
    dbtype prot."""
    url = write_source(tmp_path, PROT_FASTA)
    db_dir = tmp_path / SYN
    manifest = fetch_provider(syn_provider(url, dbtype="prot"), db_dir, fetched_at=FETCHED_AT)
    assert manifest.dbtype == "prot"
    assert (db_dir / "sequences.pin").is_file()
    assert not (db_dir / "sequences.mmi").exists()
    by_gene = {record.gene: record for record in read_records(db_dir / "records.jsonl")}
    assert by_gene["ppp_x"].sequence == SEQ_PROT_NORM


class _FakeResponse:
    """Minimal urlopen() result for the User-Agent test: read(size) hands
    out the payload, then b'' (EOF); context-manager shaped like the real
    addinfourl (the download loop uses ``with urlopen(...)``)."""

    def __init__(self, payload: bytes) -> None:
        self._remaining = payload

    def read(self, size: int = -1) -> bytes:
        take = len(self._remaining) if size < 0 else size
        data, self._remaining = self._remaining[:take], self._remaining[take:]
        return data

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_exc: object) -> None:
        return None


def test_fetch_provider_sends_gapit_user_agent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a monkeypatched urlopen standing in for the network, When a
    provider is fetched, Then the outbound call receives a urllib Request
    whose User-Agent identifies gapit (mgc.ac.cn 403s 'Python-urllib' UAs
    specifically — Wave E diagnosis) and the pipeline still completes from
    the streamed body."""
    captured: dict[str, str] = {}

    def fake_urlopen(request: urllib.request.Request) -> _FakeResponse:
        captured["url"] = request.full_url
        captured["ua"] = request.get_header("User-agent") or ""
        return _FakeResponse(UPSTREAM_FASTA.encode("utf-8"))

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    url = write_source(tmp_path)
    manifest = fetch_provider(syn_provider(url), tmp_path / SYN, fetched_at=FETCHED_AT)

    assert manifest.n_records == 3  # the fake body fed the whole pipeline
    assert captured["url"] == url
    assert "gapit/" in captured["ua"]
