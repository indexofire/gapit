"""Tests for bundled database snapshots (Wave G): make_snapshot determinism,
snapshot-path install via fetch_provider, and from_source routing.

Offline by construction: snapshot archives are BUILT IN-TEST via
``make_snapshot`` (no committed binaries — card/vfdb precedent), and the
bundled-archive lookup is exercised through the
``gapit.providers.common._snapshot_path`` seam (monkeypatched via string
setattr to a tmp archive) — never by monkeypatching importlib itself.
"""

import io
import tarfile
from collections.abc import Iterable
from pathlib import Path

import pytest

from gapit.errors import DatabaseError
from gapit.fasta import iter_fasta
from gapit.providers.common import Provider, fetch_provider
from gapit.providers.snapshots import make_snapshot
from gapit.records import (
    Manifest,
    Record,
    read_records,
    write_manifest,
    write_records,
)

SYN = "synamr"
SNAPSHOT_FILE = "synamr.tar.gz"
FETCHED_AT = "2026-09-18T09:00:00Z"
FRESH_AT = "2026-09-18T10:00:00Z"
SEQ_A = "ACGTAG" * 10  # 60 bp: minimap2-indexable (Wave A3 pin)
SEQ_B = "ACGGCTAGAT" * 6
UPSTREAM_FASTA = f">syn_a first synthetic gene\n{SEQ_A}\n>syn_b second synthetic gene\n{SEQ_B}\n"


def syn_transform(workdir: Path) -> Iterable[Record]:
    """Load the downloaded upstream fasta and yield Records (db = SYN)."""
    for fasta in iter_fasta(workdir / "upstream.fa"):
        yield Record(db=SYN, gene=fasta.id, sequence=fasta.sequence)


def syn_provider(*, snapshot: str | None, url: str) -> Provider:
    """The synthetic provider, optionally declaring a bundled snapshot."""
    return Provider(
        name=SYN,
        description="synthetic Wave G test provider",
        source_urls=(url,),
        dbtype="nucl",
        transform=syn_transform,
        snapshot=snapshot,
    )


def write_snapshot_source(dest: Path, records: tuple[Record, ...], upstream_version: str) -> None:
    """Hand-build a minimal INSTALLED db (records.jsonl + gapit-manifest.json)
    — the directory shape make_snapshot archives."""
    dest.mkdir(parents=True, exist_ok=True)
    write_records(records, dest / "records.jsonl")
    write_manifest(
        Manifest(
            name=dest.name,
            source_urls=(),
            fetched_at="2020-01-01T00:00:00Z",
            sha256="0" * 64,
            n_records=len(records),
            dbtype="nucl",
            upstream_version=upstream_version,
        ),
        dest / "gapit-manifest.json",
    )


def build_snapshot(tmp_path: Path, records: tuple[Record, ...], upstream_version: str) -> Path:
    """A valid snapshot archive in tmp_path built from a synthetic db dir."""
    source = tmp_path / "snapshotted"
    write_snapshot_source(source, records, upstream_version)
    archive = tmp_path / SNAPSHOT_FILE
    make_snapshot(source, archive)
    return archive


def patch_snapshot(monkeypatch: pytest.MonkeyPatch, archive: Path) -> None:
    """Point the common._snapshot_path seam at ``archive`` for any provider
    that declares a snapshot filename (None otherwise)."""

    def fake_snapshot_path(provider: Provider) -> Path | None:
        return archive if provider.snapshot is not None else None

    monkeypatch.setattr("gapit.providers.common._snapshot_path", fake_snapshot_path)


def test_make_snapshot_is_deterministic(tmp_path: Path) -> None:
    """Given an installed synthetic db, When snapshotted twice, Then the two
    archives are byte-identical (sorted names, PAX, gzip mtime 0) and contain
    exactly records.jsonl + gapit-manifest.json, in sorted order, whose bytes
    equal the source files."""
    source = tmp_path / "installed"
    write_snapshot_source(
        source, (Record(db=SYN, gene="syn_a", sequence=SEQ_A),), "syn-1.0"
    )

    first = tmp_path / "first.tar.gz"
    second = tmp_path / "second.tar.gz"
    make_snapshot(source, first)
    make_snapshot(source, second)

    assert first.read_bytes() == second.read_bytes()
    with tarfile.open(first, "r:gz") as tar:
        assert tar.getnames() == ["gapit-manifest.json", "records.jsonl"]
        extracted = tar.extractfile("records.jsonl")
        assert extracted is not None
        assert extracted.read() == (source / "records.jsonl").read_bytes()


def test_snapshot_install_matches_records_direct_build(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """Given a valid snapshot archive resolved via the _snapshot_path seam,
    When fetched without touching the network, Then records.jsonl and
    sequences bytes equal a records-direct network build over the same
    records, the fresh manifest inherits the ARCHIVED upstream_version while
    fetched_at is fresh, and every artifact is built."""
    source_fa = tmp_path / "upstream.fa"
    source_fa.write_text(UPSTREAM_FASTA, encoding="utf-8")
    url = source_fa.as_uri()

    reference = tmp_path / "reference"
    fetch_provider(syn_provider(snapshot=None, url=url), reference, fetched_at=FETCHED_AT)
    records = tuple(read_records(reference / "records.jsonl"))
    assert len(records) == 2

    archive = build_snapshot(tmp_path, records, "syn-1.0")
    patch_snapshot(monkeypatch, archive)

    installed = tmp_path / "installed"
    manifest = fetch_provider(
        syn_provider(snapshot=SNAPSHOT_FILE, url=url), installed, fetched_at=FRESH_AT, quiet=False
    )

    assert (installed / "records.jsonl").read_bytes() == (reference / "records.jsonl").read_bytes()
    assert (installed / "sequences").read_bytes() == (reference / "sequences").read_bytes()
    assert manifest.upstream_version == "syn-1.0"
    assert manifest.fetched_at == FRESH_AT
    assert manifest.n_records == 2
    for name in ("sequences.nin", "sequences.mmi", "gapit-manifest.json"):
        assert (installed / name).is_file(), name

    captured = capsys.readouterr()
    assert captured.out == ""  # stdout purity
    assert "snapshot" in captured.err


def test_from_source_routes_to_network(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a provider whose snapshot IS resolvable, When fetched with
    from_source=True and a dead source URL, Then DOWNLOAD_FAILED — the flag
    routed around the snapshot to the (failing) network path; a snapshot-path
    implementation would instead succeed."""
    dead_url = (tmp_path / "no-such-file.fa").as_uri()
    archive = build_snapshot(
        tmp_path, (Record(db=SYN, gene="syn_a", sequence=SEQ_A),), "syn-1.0"
    )
    patch_snapshot(monkeypatch, archive)

    with pytest.raises(DatabaseError) as excinfo:
        fetch_provider(
            syn_provider(snapshot=SNAPSHOT_FILE, url=dead_url),
            tmp_path / "db",
            fetched_at=FETCHED_AT,
            from_source=True,
        )
    error = excinfo.value
    assert error.code == "DOWNLOAD_FAILED"
    assert error.exit_code == 4
    assert error.context == {"url": dead_url}


def test_corrupt_snapshot_raises_typed_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a snapshot archive missing its gapit-manifest.json member, When
    fetched, Then DatabaseError SNAPSHOT_INVALID (no untyped KeyError escape
    — the Wave E card.py lesson)."""
    records_only = tmp_path / "records-only"
    records_only.mkdir()
    write_records(
        (Record(db=SYN, gene="syn_a", sequence=SEQ_A),), records_only / "records.jsonl"
    )
    archive = tmp_path / SNAPSHOT_FILE
    with tarfile.open(archive, "w:gz") as tar:
        info = tarfile.TarInfo("records.jsonl")
        data = (records_only / "records.jsonl").read_bytes()
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    patch_snapshot(monkeypatch, archive)

    dead_url = (tmp_path / "no-such-file.fa").as_uri()
    with pytest.raises(DatabaseError) as excinfo:
        fetch_provider(
            syn_provider(snapshot=SNAPSHOT_FILE, url=dead_url),
            tmp_path / "db",
            fetched_at=FETCHED_AT,
        )
    error = excinfo.value
    assert error.code == "SNAPSHOT_INVALID"
    assert error.exit_code == 4
    assert error.context["missing"] == "gapit-manifest.json"
