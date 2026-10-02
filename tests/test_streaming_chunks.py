"""Streaming chunk-contract tests: per-file emission through an optional
``emit`` sink, byte-identity of chunk concatenation with the buffered render,
head-of-line timing (file 1 emits before a slow file 2 completes), and the
cluster engine's twin contract.

The renderer-level byte-proofs pin the split helpers: preamble + per-file
chunks concatenated MUST equal the one-shot render, so streaming never
changes a single output byte (goldens untouched). CLI ``--output`` behavior
lives in test_cli_screen.py.
"""

import re
import shutil
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

import pytest
from typer.testing import CliRunner

from gapit import screening, screening_cluster
from gapit.blast import BlastRow
from gapit.cli import app
from gapit.cluster import ClusterParams, ClusterReport, load_features, screen_cluster_file
from gapit.db import Database, make_blast_db
from gapit.fasta import iter_fasta
from gapit.formats.cluster_md import cluster_md_file_chunk, cluster_md_head, render_cluster_md
from gapit.formats.md import md_report_chunk, md_report_preamble, render_markdown
from gapit.formats.tsv import format_tsv, tsv_file_chunk, tsv_preamble
from gapit.gbfeatures import FeaturesDocument
from gapit.hits import process_rows
from gapit.report import Report, ScreeningParams
from gapit.screening import OutputFormat, find_database, run_screen
from gapit.screening_cluster import run_cluster_screen

FIXTURE_DB_DIR = Path(__file__).parent / "data" / "db"
CONTIGS = Path(__file__).parent / "data" / "contigs"
GOLDEN = Path(__file__).parent / "golden"
MULTI_FILES = ("full.fa", "gap.fa", "none.fa", "sort.fa")
PINNED_NOW = datetime(2026, 9, 17, 12, 0, 0, tzinfo=UTC)
PARAMS = ScreeningParams(db="tinyamr")
TET_ID = "tinyamr~~~tetA~~~NC_000913.3:100-900~~~TETRACYCLINE"

runner = CliRunner()


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    """Fresh datadir with a built tinyamr index, one per test."""
    target = tmp_path / "datadir"
    shutil.copytree(FIXTURE_DB_DIR, target)
    make_blast_db(target / "tinyamr" / "sequences", "tinyamr")
    return target


def _reports() -> list[Report]:
    """Synthetic multi-file reports (one hit, one empty) — no BLAST needed."""
    row = BlastRow(
        qseqid="contig1",
        qstart=1,
        qend=79,
        qlen=79,
        sseqid=TET_ID,
        sstart=1,
        send=79,
        slen=79,
        sstrand="plus",
        evalue=1e-40,
        length=79,
        pident=100.0,
        gaps=0,
        gapopen=0,
        stitle=TET_ID + " tetracycline efflux pump TetA",
    )
    hit = process_rows([row], mincov=0.0, default_db="tinyamr")[0]
    return [
        Report(file="one.fa", hits=(hit,)),
        Report(file="two.fa", hits=()),
    ]


def _run(
    datadir: Path,
    output_format: OutputFormat,
    files: list[Path] | None = None,
    jobs: int = 1,
) -> str:
    return run_screen(
        files or [CONTIGS / name for name in MULTI_FILES],
        "tinyamr",
        datadir,
        80.0,
        80.0,
        1,
        jobs,
        None,
        quiet=True,
        noheader=False,
        nopath=True,
        debug=False,
        output_format=output_format,
    )


# ------------------------------------------------- renderer byte-proofs --


@pytest.mark.parametrize(("noheader", "csv"), [(False, False), (True, False), (False, True)])
def test_tsv_chunks_concatenate_to_format_tsv(noheader: bool, csv: bool) -> None:
    """Given reports, When split into preamble + per-file chunks, Then the
    concatenation is byte-identical to format_tsv's one-shot render."""
    reports = _reports()
    chunks = [tsv_preamble(csv=csv, noheader=noheader)]
    chunks += [tsv_file_chunk(report, csv=csv, nopath=False) for report in reports]
    assert "".join(chunks) == format_tsv(reports, csv=csv, noheader=noheader, nopath=False)


def test_md_chunks_concatenate_to_render_markdown() -> None:
    """Given reports, When split into preamble + per-file sections, Then the
    concatenation is byte-identical to render_markdown's one-shot render
    (the frontmatter is static, so the streamer emits it first)."""
    reports = _reports()
    preamble = md_report_preamble(PARAMS, now=PINNED_NOW)
    assert preamble + "".join(md_report_chunk(report) for report in reports) == render_markdown(
        reports, PARAMS, now=PINNED_NOW
    )


# ------------------------------------------- use-case emit contract --


@pytest.mark.parametrize(
    ("output_format", "golden_name"),
    [
        (OutputFormat.tsv, "tinyamr_multi_nopath.tsv"),
        (OutputFormat.csv, "tinyamr_multi_nopath.csv"),
    ],
)
def test_emit_chunks_match_buffered_render_and_golden(
    datadir: Path, output_format: OutputFormat, golden_name: str
) -> None:
    """Given a streamed run (emit=collector) beside a buffered run, When
    compared, Then chunk concatenation, the returned string, and the
    committed golden are all byte-identical — streaming changes nothing."""
    collected: list[str] = []
    streamed = run_screen(
        [CONTIGS / name for name in MULTI_FILES],
        "tinyamr",
        datadir,
        80.0,
        80.0,
        1,
        1,
        None,
        quiet=True,
        noheader=False,
        nopath=True,
        debug=False,
        output_format=output_format,
        emit=collected.append,
    )
    buffered = _run(datadir, output_format)
    assert streamed == buffered
    assert "".join(collected) == buffered
    assert buffered == (GOLDEN / golden_name).read_text(encoding="utf-8")


def _without_timestamps(rendered: str) -> str:
    """created_at differs across runs on a second boundary; comparisons
    normalize it away (the quiet-comparison precedent)."""
    return re.sub(r"created_at: \S+|\"created_at\": \"[^\"]+\"", "<ts>", rendered)


def test_emit_json_single_document_matches_buffered(datadir: Path) -> None:
    """Given json (the single-document contract), When streamed through
    emit, Then the one chunk equals the buffered render."""
    collected: list[str] = []
    streamed = run_screen(
        [CONTIGS / name for name in MULTI_FILES],
        "tinyamr",
        datadir,
        80.0,
        80.0,
        1,
        1,
        None,
        quiet=True,
        noheader=False,
        nopath=True,
        debug=False,
        output_format=OutputFormat.json,
        emit=collected.append,
    )
    buffered = _run(datadir, OutputFormat.json)
    assert _without_timestamps(streamed) == _without_timestamps(buffered)
    assert _without_timestamps("".join(collected)) == _without_timestamps(buffered)
    assert len(collected) == 1


def test_emit_tsv_chunk_arrives_per_file_and_in_order(datadir: Path) -> None:
    """Given a 4-file tsv run with emit, When the chunk sequence is
    inspected, Then the header lands first and each file's rows arrive as
    their own chunk in input order (one chunk per file after the preamble)."""
    collected: list[str] = []
    run_screen(
        [CONTIGS / name for name in MULTI_FILES],
        "tinyamr",
        datadir,
        80.0,
        80.0,
        1,
        1,
        None,
        quiet=True,
        noheader=False,
        nopath=True,
        debug=False,
        output_format=OutputFormat.tsv,
        emit=collected.append,
    )
    assert collected[0].startswith("#FILE\t")
    # none.fa has zero hits, so it emits no chunk (empty chunks are skipped);
    # every file WITH hits arrives as its own chunk, in input order.
    assert [chunk.splitlines()[0].split("\t")[0] for chunk in collected[1:]] == [
        "full.fa",
        "gap.fa",
        "sort.fa",
    ]


def test_file1_emits_before_slow_file2_completes(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given --jobs 2 with file 2 artificially slowed, When streaming, Then
    file 1's rows are emitted while file 2 is still screening — the
    head-of-line proof (emit fires per completed file, not at run end)."""
    real = screening.screen_file
    gap_done: dict[str, float] = {}

    def slow_then_record(
        query: Path,
        database: Database,
        params: ScreeningParams,
        *,
        dbtype: Literal["nucl", "prot"],
        debug: bool = False,
        merge_fragments: bool = False,
    ) -> Report:
        report = real(
            query, database, params, dbtype=dbtype, debug=debug, merge_fragments=merge_fragments
        )
        if query.name == "gap.fa":
            time.sleep(0.4)
            gap_done["at"] = time.monotonic()
        return report

    monkeypatch.setattr(screening, "screen_file", slow_then_record)
    stamps: list[tuple[str, float]] = []
    run_screen(
        [CONTIGS / "full.fa", CONTIGS / "gap.fa"],
        "tinyamr",
        datadir,
        80.0,
        80.0,
        1,
        2,
        None,
        quiet=True,
        noheader=False,
        nopath=True,
        debug=False,
        output_format=OutputFormat.tsv,
        emit=lambda chunk: stamps.append((chunk, time.monotonic())),
    )
    full_chunk, full_at = next((chunk, at) for chunk, at in stamps if chunk.startswith("full.fa\t"))
    assert full_chunk.splitlines()[0].split("\t")[0] == "full.fa"
    assert full_at < gap_done["at"]


# ------------------------------------------------------- md streaming --


def test_emit_md_chunk_arrives_per_file_and_in_order(datadir: Path) -> None:
    """Given a 4-file md run with emit, When the chunk sequence is
    inspected, Then the static frontmatter lands first and every file's
    section (zero-hit files included) arrives as its own chunk in input
    order — md streams exactly like the tsv chunks."""
    collected: list[str] = []
    run_screen(
        [CONTIGS / name for name in MULTI_FILES],
        "tinyamr",
        datadir,
        80.0,
        80.0,
        1,
        1,
        None,
        quiet=True,
        noheader=False,
        nopath=True,
        debug=False,
        output_format=OutputFormat.md,
        emit=collected.append,
    )
    assert collected[0].startswith("---\nschema: gapit.report/1\n")
    assert "# gapit screening report" in collected[0]
    # md headings keep the path as given (--nopath rewrites only the TSV
    # FILE column), so the expected first lines carry the fixture paths.
    assert [chunk.splitlines()[0] for chunk in collected[1:]] == [
        f"## `{CONTIGS / name}`" for name in MULTI_FILES
    ]
    assert "files:" not in collected[0] and "hits:" not in collected[0]


@pytest.mark.parametrize("output_format", [OutputFormat.tsv, OutputFormat.md])
def test_file1_emits_before_slow_file2_completes_across_formats(
    datadir: Path, monkeypatch: pytest.MonkeyPatch, output_format: OutputFormat
) -> None:
    """Given --jobs 2 with file 2 artificially slowed, When streaming tsv
    or md, Then file 1's chunk is emitted while file 2 is still screening —
    the md frontmatter is static, so nothing waits for the last file."""
    real = screening.screen_file
    gap_done: dict[str, float] = {}

    def slow_then_record(
        query: Path,
        database: Database,
        params: ScreeningParams,
        *,
        dbtype: Literal["nucl", "prot"],
        debug: bool = False,
        merge_fragments: bool = False,
    ) -> Report:
        report = real(
            query, database, params, dbtype=dbtype, debug=debug, merge_fragments=merge_fragments
        )
        if query.name == "gap.fa":
            time.sleep(0.4)
            gap_done["at"] = time.monotonic()
        return report

    monkeypatch.setattr(screening, "screen_file", slow_then_record)
    stamps: list[tuple[str, float]] = []
    run_screen(
        [CONTIGS / "full.fa", CONTIGS / "gap.fa"],
        "tinyamr",
        datadir,
        80.0,
        80.0,
        1,
        2,
        None,
        quiet=True,
        noheader=False,
        nopath=True,
        debug=False,
        output_format=output_format,
        emit=lambda chunk: stamps.append((chunk, time.monotonic())),
    )
    file1_heading = f"## `{CONTIGS / 'full.fa'}`"
    file1_at = next(
        at
        for chunk, at in stamps
        # tsv chunks start "full.fa\t" (nopath); md headings keep the path.
        if chunk.startswith("full.fa\t") or chunk.splitlines()[0] == file1_heading
    )
    assert file1_at < gap_done["at"]


# ----------------------------------------------------- cluster engine --


def _cluster_samples(datadir: Path, tmp_path: Path) -> list[Path]:
    """Build the cluster fixture db and two derived sample assemblies."""
    data = Path(__file__).parent / "data" / "cluster"
    built = runner.invoke(
        app, ["db", "build", "cps", str(data / "screening.gbk"), "--datadir", str(datadir)]
    )
    assert built.exit_code == 0, built.stderr
    samples = tmp_path / "samples"
    samples.mkdir()
    locus_a = {record.id: record.sequence for record in iter_fasta(datadir / "cps" / "sequences")}[
        "locusA"
    ]
    (samples / "a.fa").write_text(f">ctg_a\n{locus_a}\n", encoding="utf-8")
    (samples / "b.fa").write_text(f">ctg_b\n{locus_a[:800]}\n", encoding="utf-8")
    return [samples / "a.fa", samples / "b.fa"]


def test_cluster_emit_chunks_match_buffered_render(tmp_path: Path) -> None:
    """Given the cluster fixture db and two samples, When streamed through
    emit beside a buffered run, Then chunk concatenation equals the buffered
    render byte-for-byte and each file's row arrives as its own chunk."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    files = _cluster_samples(datadir, tmp_path)
    database = find_database(datadir, "cps")
    collected: list[str] = []
    streamed = run_cluster_screen(
        files,
        database,
        ClusterParams(db="cps"),
        output_format=OutputFormat.tsv,
        noheader=False,
        nopath=True,
        quiet=True,
        emit=collected.append,
    )
    buffered = run_cluster_screen(
        files,
        database,
        ClusterParams(db="cps"),
        output_format=OutputFormat.tsv,
        noheader=False,
        nopath=True,
        quiet=True,
    )
    assert streamed == buffered
    assert "".join(collected) == buffered
    assert [chunk.splitlines()[0].split("\t")[0] for chunk in collected[1:]] == ["a.fa", "b.fa"]


def test_cluster_md_head_and_chunks_concatenate_to_render(tmp_path: Path) -> None:
    """Given cluster reports, When split into static head + per-file
    sections, Then the concatenation is byte-identical to render_cluster_md
    (the summary row rides inside each file's own section)."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    files = _cluster_samples(datadir, tmp_path)
    database = find_database(datadir, "cps")
    params = ClusterParams(db="cps")
    reports = [
        screen_cluster_file(path, database, load_features(database), params) for path in files
    ]
    assert cluster_md_head(params, now=PINNED_NOW) + "".join(
        cluster_md_file_chunk(report, typed=False) for report in reports
    ) == render_cluster_md(reports, params, now=PINNED_NOW, typed=False)


def test_cluster_md_emits_per_file_with_own_summary_row(tmp_path: Path) -> None:
    """Given the cluster fixture db and two samples streamed as md, When
    the chunk sequence is inspected, Then the static frontmatter lands
    first and each file's section (heading + its OWN one-row summary table
    + gene table) arrives as its own chunk — no all-files summary, no
    end-buffering."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    files = _cluster_samples(datadir, tmp_path)
    database = find_database(datadir, "cps")
    collected: list[str] = []
    streamed = run_cluster_screen(
        files,
        database,
        ClusterParams(db="cps"),
        output_format=OutputFormat.md,
        noheader=False,
        nopath=True,
        quiet=True,
        emit=collected.append,
    )
    buffered = run_cluster_screen(
        files,
        database,
        ClusterParams(db="cps"),
        output_format=OutputFormat.md,
        noheader=False,
        nopath=True,
        quiet=True,
    )
    assert _without_timestamps(streamed) == _without_timestamps(buffered)
    assert _without_timestamps("".join(collected)) == _without_timestamps(buffered)
    assert collected[0].startswith("---\nschema: gapit.cluster/1\n")
    assert "files:" not in collected[0]
    assert [chunk.splitlines()[0] for chunk in collected[1:]] == [f"## `{file}`" for file in files]
    for chunk, file in zip(collected[1:], files, strict=True):
        assert f"| {file} |" in chunk  # the summary row names its file


def test_cluster_md_file1_emits_before_slow_file2_completes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a cluster md run with file 2 artificially slowed, When
    streaming, Then file 1's section is emitted while file 2 is still
    screening — the md path never waits for the last file."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    files = _cluster_samples(datadir, tmp_path)
    database = find_database(datadir, "cps")
    real = screening_cluster.screen_cluster_file
    b_done: dict[str, float] = {}

    def slow_b(
        path: Path,
        database: Database,
        features: FeaturesDocument,
        params: ClusterParams,
        *,
        debug: bool = False,
    ) -> ClusterReport:
        report = real(path, database, features, params, debug=debug)
        if path.name == "b.fa":
            time.sleep(0.4)
            b_done["at"] = time.monotonic()
        return report

    monkeypatch.setattr(screening_cluster, "screen_cluster_file", slow_b)
    stamps: list[tuple[str, float]] = []
    run_cluster_screen(
        files,
        database,
        ClusterParams(db="cps"),
        output_format=OutputFormat.md,
        noheader=False,
        nopath=True,
        quiet=True,
        emit=lambda chunk: stamps.append((chunk, time.monotonic())),
    )
    a_heading = f"## `{files[0]}`"
    a_chunk, a_at = next((chunk, at) for chunk, at in stamps if chunk.splitlines()[0] == a_heading)
    assert f"| {files[0]} |" in a_chunk
    assert a_at < b_done["at"]
