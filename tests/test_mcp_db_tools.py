"""MCP db tool tests: db_fetch / db_build / db_search / db_outdated driven
through the real serve loop.

Same offline patterns as the CLI suites: hand-written manifest + records
databases for the read-only queries (test_cli_db_query.py), file://
providers and in-test snapshot archives for fetch (test_cli_db_fetch.py),
and real makeblastdb/blastn for the build-then-screen agent story
(test_cli_db_build.py). No network; tmp_path datadirs only.
"""

import json
import random
from collections.abc import Iterable
from datetime import UTC, datetime, timedelta
from io import StringIO
from pathlib import Path
from typing import Any

import pytest

from gapit.fasta import iter_fasta
from gapit.mcp import serve
from gapit.providers.common import Provider
from gapit.providers.snapshots import make_snapshot
from gapit.records import Manifest, Record, write_manifest, write_records

SYN = "synamr"
SEQ = "ACGTAG" * 10  # 60 bp, the synthetic-fixture convention
SEQ_B = "ACGGCTAGAT" * 6  # distinct second sequence (dedupe would drop a twin)
UPSTREAM_FASTA = f">syn_a first synthetic gene\n{SEQ}\n>syn_b second synthetic gene\n{SEQ_B}\n"
# Distinct seeds share no 11-mer in practice, so the build fixture's genes
# cannot cross-match under blastn (test_cli_db_build.py note).
SEQ_A = "".join(random.Random(1).choice("ACGT") for _ in range(240))


def exchange(*lines_or_messages: object) -> list[dict[str, Any]]:
    """Feed the loop one line per item; return the parsed response objects."""
    raw = "".join(
        (item if isinstance(item, str) else json.dumps(item)) + "\n" for item in lines_or_messages
    )
    out = StringIO()
    serve(StringIO(raw), out)
    return [json.loads(line) for line in out.getvalue().splitlines()]


def tool_call(name: str, arguments: dict[str, object], *, id_value: int = 1) -> dict[str, object]:
    return {
        "jsonrpc": "2.0",
        "id": id_value,
        "method": "tools/call",
        "params": {"name": name, "arguments": arguments},
    }


def call_tool(name: str, arguments: dict[str, object]) -> tuple[bool, str]:
    """One tools/call round trip -> (isError, text)."""
    (response,) = exchange(tool_call(name, arguments))
    result = response["result"]
    return result["isError"], result["content"][0]["text"]


def envelope_code(text: str) -> str:
    return json.loads(text)["code"]


def make_db(
    datadir: Path,
    name: str,
    fetched_at: str,
    records: tuple[Record, ...],
    *,
    upstream_version: str = "",
) -> None:
    """Hand-write an installed database: manifest + records.jsonl."""
    db_dir = datadir / name
    db_dir.mkdir(parents=True)
    write_records(records, db_dir / "records.jsonl")
    write_manifest(
        Manifest(
            name=name,
            source_urls=(),
            fetched_at=fetched_at,
            sha256="0" * 64,
            n_records=len(records),
            dbtype="nucl",
            upstream_version=upstream_version,
        ),
        db_dir / "gapit-manifest.json",
    )


def rec(
    db: str,
    gene: str,
    *,
    accession: str = "",
    function: tuple[str, ...] = (),
    product: str = "n/a",
) -> Record:
    """One synthetic record (60 bp sequence)."""
    return Record(
        db=db, gene=gene, sequence=SEQ, accession=accession, function=function, product=product
    )


def query_datadir(tmp_path: Path) -> Path:
    """Datadir with two searchable databases (3 'tet' genes across them)."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fresh_at = (datetime.now(UTC) - timedelta(days=1)).strftime("%Y-%m-%dT%H:%M:%SZ")
    make_db(
        datadir,
        "atiny",
        fresh_at,
        (
            rec("atiny", "tetA", accession="ACC-1", function=("TETRACYCLINE",), product="efflux"),
            rec("atiny", "tetM", accession="ACC-2", product="ribosomal protection"),
        ),
    )
    make_db(
        datadir, "btiny", "2020-01-01T00:00:00Z", (rec("btiny", "tetX", product="monoxygenase"),)
    )
    return datadir


# ------------------------------------------------------------- tools/list --


def test_tools_list_carries_the_four_db_tools_with_schemas() -> None:
    """Given tools/list, When served, Then db_fetch/db_build/db_search/
    db_outdated are advertised with object schemas and correct requireds."""
    (response,) = exchange({"jsonrpc": "2.0", "id": 1, "method": "tools/list"})
    by_name = {tool["name"]: tool for tool in response["result"]["tools"]}
    assert by_name["db_fetch"]["inputSchema"]["required"] == []
    assert by_name["db_build"]["inputSchema"]["required"] == ["name", "fasta"]
    assert by_name["db_search"]["inputSchema"]["required"] == ["term"]
    assert by_name["db_outdated"]["inputSchema"]["required"] == []
    for name in ("db_fetch", "db_build", "db_search", "db_outdated"):
        assert by_name[name]["description"]
        assert by_name[name]["inputSchema"]["type"] == "object"


# ---------------------------------------------------------------- search --


def test_db_search_returns_tsv_hits_sorted_and_truncated(tmp_path: Path) -> None:
    """Given a two-db datadir with three 'tet' genes, When db_search('tet')
    with limit 2, Then exactly two DB-sorted TSV rows come back (columns
    DB/GENE/ACCESSION/FUNCTION/PRODUCT/LENGTH); without limit all three do."""
    datadir = query_datadir(tmp_path)
    is_error, text = call_tool("db_search", {"term": "tet", "limit": 2, "datadir": str(datadir)})
    assert is_error is False
    lines = text.splitlines()
    assert len(lines) == 2
    assert lines[0].startswith("atiny\ttetA\tACC-1\tTETRACYCLINE\tefflux\t60")
    is_error, text = call_tool("db_search", {"term": "tet", "datadir": str(datadir)})
    assert is_error is False
    assert [line.split("\t")[1] for line in text.splitlines()] == ["tetA", "tetM", "tetX"]


def test_db_search_field_exact_and_unknown_db(tmp_path: Path) -> None:
    """Given the fixture datadir, When db_search uses field=gene + exact,
    Then only the exact gene matches; an unknown db and a bad field are
    USAGE_ERROR envelopes via isError."""
    datadir = query_datadir(tmp_path)
    is_error, text = call_tool(
        "db_search", {"term": "tetA", "field": "gene", "exact": True, "datadir": str(datadir)}
    )
    assert is_error is False
    assert [line.split("\t")[1] for line in text.splitlines()] == ["tetA"]
    is_error, text = call_tool("db_search", {"term": "tet", "db": "nope", "datadir": str(datadir)})
    assert is_error is True
    assert envelope_code(text) == "USAGE_ERROR"
    is_error, text = call_tool("db_search", {"term": "tet", "field": "bogus"})
    assert is_error is True
    assert envelope_code(text) == "USAGE_ERROR"


# -------------------------------------------------------------- outdated --


def test_db_outdated_reports_tsv_table(tmp_path: Path) -> None:
    """Given one fresh and one stale database, When db_outdated, Then the
    NAME/FETCHED_AT/AGE_DAYS/STATUS table text comes back (fresh=ok,
    stale=stale); days=0 marks everything stale."""
    datadir = query_datadir(tmp_path)
    is_error, text = call_tool("db_outdated", {"datadir": str(datadir)})
    assert is_error is False
    lines = text.splitlines()
    assert lines[0] == "NAME\tFETCHED_AT\tAGE_DAYS\tSTATUS"
    table = {line.split("\t")[0]: line.split("\t") for line in lines[1:]}
    assert table["atiny"][3] == "ok"
    assert table["btiny"][3] == "stale"
    is_error, text = call_tool("db_outdated", {"days": 0, "datadir": str(datadir)})
    assert is_error is False
    assert all(line.split("\t")[3] == "stale" for line in text.splitlines()[1:])


def test_db_outdated_empty_datadir_is_database_error(tmp_path: Path) -> None:
    """Given a datadir with no installed databases, When db_outdated, Then
    isError with the DATADIR_EMPTY envelope (CLI exit-4 semantics)."""
    datadir = tmp_path / "empty"
    datadir.mkdir()
    is_error, text = call_tool("db_outdated", {"datadir": str(datadir)})
    assert is_error is True
    assert envelope_code(text) == "DATADIR_EMPTY"


# ----------------------------------------------------------------- build --


def test_db_build_then_screen_full_agent_story(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a local FASTA and a fresh datadir on $GAPIT_DATADIR, When
    db_build, Then a JSON receipt returns and the native artifacts exist;
    and When screen runs against the new database in the SAME session, Then
    the gene is found — the self-provisioning agent story."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    monkeypatch.setenv("GAPIT_DATADIR", str(datadir))
    fasta = tmp_path / "my_genes.fa"
    fasta.write_text(f">demov2 demo beta-lactamase variant 2\n{SEQ_A}\n", encoding="utf-8")
    query = tmp_path / "query.fa"
    query.write_text(f">contig1\n{SEQ_A[:200]}\n", encoding="utf-8")

    is_error, text = call_tool("db_build", {"name": "myamr", "fasta": str(fasta)})
    assert is_error is False
    assert json.loads(text) == {
        "db": "myamr",
        "records": 1,
        "dbtype": "nucl",
        "destination": str(datadir / "myamr"),
    }
    for artifact in ("sequences", "sequences.nin", "gapit-manifest.json", "records.jsonl"):
        assert (datadir / "myamr" / artifact).is_file(), artifact

    is_error, text = call_tool("screen", {"files": [str(query)], "db": "myamr"})
    assert is_error is False
    document = json.loads(text)
    assert document["schema"] == "gapit.report/1"
    (hit,) = document["files"][0]["hits"]
    assert hit["gene"] == "demov2"


def test_db_build_refuses_overwrite_without_force(tmp_path: Path) -> None:
    """Given an already-built database, When db_build runs again without
    force, Then isError DB_ALREADY_EXISTS; with force=true it rebuilds."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fasta = tmp_path / "my_genes.fa"
    fasta.write_text(f">demov2 demo gene\n{SEQ_A}\n", encoding="utf-8")
    args: dict[str, object] = {"name": "myamr", "fasta": str(fasta), "datadir": str(datadir)}
    assert call_tool("db_build", args)[0] is False
    is_error, text = call_tool("db_build", args)
    assert is_error is True
    assert envelope_code(text) == "DB_ALREADY_EXISTS"
    is_error, text = call_tool("db_build", {**args, "force": True})
    assert is_error is False
    assert json.loads(text)["records"] == 1


def test_db_build_missing_fasta_is_input_error(tmp_path: Path) -> None:
    """Given a FASTA path that does not exist, When db_build, Then isError
    INPUT_NOT_FOUND naming the file (CLI semantics)."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fasta = tmp_path / "absent.fa"
    is_error, text = call_tool(
        "db_build", {"name": "myamr", "fasta": str(fasta), "datadir": str(datadir)}
    )
    assert is_error is True
    assert envelope_code(text) == "INPUT_NOT_FOUND"


# ----------------------------------------------------------------- fetch --


def syn_transform(workdir: Path) -> Iterable[Record]:
    for fasta in iter_fasta(workdir / "upstream.fa"):
        yield Record(db=SYN, gene=fasta.id, sequence=fasta.sequence)


def patch_file_registry(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """One file:// provider in gapit.cmd_db.REGISTRY + an empty datadir."""
    source = tmp_path / "upstream.fa"
    source.write_text(UPSTREAM_FASTA, encoding="utf-8")
    monkeypatch.setattr(
        "gapit.cmd_db.REGISTRY",
        {
            SYN: Provider(
                name=SYN,
                description="synthetic provider for the MCP db_fetch tests",
                source_urls=(source.as_uri(),),
                dbtype="nucl",
                transform=syn_transform,
            )
        },
    )
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    return datadir


def test_db_fetch_returns_receipt_and_builds_artifacts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a file:// provider, When db_fetch(NAME), Then one JSON receipt
    line comes back naming the destination and the database artifacts exist
    under the datadir argument."""
    datadir = patch_file_registry(tmp_path, monkeypatch)
    is_error, text = call_tool("db_fetch", {"name": SYN, "datadir": str(datadir)})
    assert is_error is False
    assert json.loads(text) == {
        "db": SYN,
        "records": 2,
        "dbtype": "nucl",
        "destination": str(datadir / SYN),
    }
    for artifact in ("sequences", "sequences.nin", "gapit-manifest.json"):
        assert (datadir / SYN / artifact).is_file(), artifact


def test_db_fetch_exists_then_force(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given an already-fetched database, When db_fetch runs again, Then
    isError DB_ALREADY_EXISTS; force=true rebuilds with the same receipt."""
    datadir = patch_file_registry(tmp_path, monkeypatch)
    args: dict[str, object] = {"name": SYN, "datadir": str(datadir)}
    assert call_tool("db_fetch", args)[0] is False
    is_error, text = call_tool("db_fetch", args)
    assert is_error is True
    assert envelope_code(text) == "DB_ALREADY_EXISTS"
    is_error, text = call_tool("db_fetch", {**args, "force": True})
    assert is_error is False
    assert json.loads(text)["records"] == 2


def test_db_fetch_unknown_name_is_usage_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a patched registry, When db_fetch names a provider not in it,
    Then isError with the USAGE_ERROR envelope (CLI exit-2 semantics)."""
    datadir = patch_file_registry(tmp_path, monkeypatch)
    is_error, text = call_tool("db_fetch", {"name": "ncbi", "datadir": str(datadir)})
    assert is_error is True
    assert envelope_code(text) == "USAGE_ERROR"
    assert SYN in text


def test_db_fetch_default_set_installs_card_vfdb_from_snapshots(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given card+vfdb providers whose snapshots carry one record each and
    whose network URLs are dead, When db_fetch with no name, Then both
    defaults install from the snapshots: two receipt lines, card first."""
    dead = (tmp_path / "dead.fa").as_uri()
    archives: dict[str, Path] = {}
    registry: dict[str, Provider] = {}
    for name in ("card", "vfdb"):
        source = tmp_path / f"{name}-source"
        source.mkdir()
        write_records(
            (Record(db=name, gene=f"snap_{name}_gene", sequence=SEQ),),
            source / "records.jsonl",
        )
        write_manifest(
            Manifest(
                name=name,
                source_urls=(),
                fetched_at="2020-01-01T00:00:00Z",
                sha256="0" * 64,
                n_records=1,
                dbtype="nucl",
                upstream_version=f"{name}-4.0",
            ),
            source / "gapit-manifest.json",
        )
        archive = tmp_path / f"{name}.tar.gz"
        make_snapshot(source, archive)
        archives[archive.name] = archive
        registry[name] = Provider(
            name=name,
            description=f"synthetic snapshot-backed {name} provider",
            source_urls=(dead,),
            dbtype="nucl",
            transform=syn_transform,
            snapshot=f"{name}.tar.gz",
        )
    monkeypatch.setattr("gapit.cmd_db.REGISTRY", registry)

    def fake_snapshot_path(provider: Provider) -> Path | None:
        return archives.get(provider.snapshot) if provider.snapshot is not None else None

    monkeypatch.setattr("gapit.providers.common._snapshot_path", fake_snapshot_path)
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    is_error, text = call_tool("db_fetch", {"datadir": str(datadir)})
    assert is_error is False
    receipts = [json.loads(line) for line in text.splitlines()]
    assert [receipt["db"] for receipt in receipts] == ["card", "vfdb"]
    assert all(receipt["records"] == 1 for receipt in receipts)
    for name in ("card", "vfdb"):
        assert (datadir / name / "gapit-manifest.json").is_file(), name
