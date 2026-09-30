"""Cluster-engine suites (stage 2): CLI dispatch, guards, goldens, MCP, and
the real-data smoke.

The fixture database is built per test from tests/data/cluster/screening.gbk
(repo convention: commit the source GBK only, build into tmp_path). Samples
are derived in-test from the built locus FASTA, so goldens stay stable:
locusA (3 genes incl. two minus-strand) exact / fragmented 60-40 across two
contigs / ~5% mutated / truncated 2-of-3 genes, plus an unrelated contig.
Real minimap2 + makeblastdb from the pixi env (test_cli_reads.py recipe).
"""

import json
import shutil
from datetime import UTC, datetime
from io import StringIO
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner, Result

from gapit.cli import app
from gapit.cluster import ClusterParams, load_features, screen_cluster_file
from gapit.db import discover_databases
from gapit.fasta import iter_fasta
from gapit.formats.cluster import format_cluster_tsv, render_cluster_json
from gapit.formats.cluster_md import render_cluster_md
from gapit.mcp import serve

DATA = Path(__file__).parent / "data" / "cluster"
GENE_DB = Path(__file__).parent / "data" / "db"
GOLDEN = Path(__file__).parent / "golden"
DB = "cps"
PINNED_NOW = datetime(2026, 9, 30, 12, 0, 0, tzinfo=UTC)
PARAMS = ClusterParams(db=DB)
REAL_CPS = Path("/tmp/opencode/vpautils/src/vpautils/database/serotype/CPSgc.gbk")

runner = CliRunner()
BASES = "ACGT"


def lcg_seq(length: int, seed: int) -> str:
    """Deterministic synthetic DNA (the fixture generator's rule)."""
    state = seed
    out: list[str] = []
    for _ in range(length):
        state = (state * 1103515245 + 12345) % (1 << 31)
        out.append(BASES[(state >> 16) % 4])
    return "".join(out)


@pytest.fixture()
def cluster_db(tmp_path: Path) -> Path:
    """Build the screening fixture db into a fresh datadir and return it."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    result = runner.invoke(
        app, ["db", "build", DB, str(DATA / "screening.gbk"), "--datadir", str(datadir)]
    )
    assert result.exit_code == 0, result.stderr
    return datadir


def _loci(datadir: Path) -> dict[str, str]:
    return {record.id: record.sequence for record in iter_fasta(datadir / DB / "sequences")}


def _write_samples(datadir: Path, tmp_path: Path) -> dict[str, Path]:
    """Derive the five sample assemblies from the built locus sequences."""
    locus_a = _loci(datadir)["locusA"]
    samples = tmp_path / "samples"
    samples.mkdir()
    mutated = "".join(
        BASES[(BASES.index(base) + 1) % 4] if offset % 20 == 0 else base
        for offset, base in enumerate(locus_a)
    )
    contents = {
        "exact.fa": f">ctg_exact\n{locus_a}\n",
        "fragmented.fa": f">ctg_left\n{locus_a[:720]}\n>ctg_right\n{locus_a[720:]}\n",
        "partial.fa": f">ctg_partial\n{locus_a[:800]}\n",
        "mutated.fa": f">ctg_mutated\n{mutated}\n",
        "unrelated.fa": f">ctg_unrelated\n{lcg_seq(1100, 999)}\n",
    }
    for name, text in contents.items():
        (samples / name).write_text(text, encoding="utf-8")
    return {path.stem: path for path in sorted(samples.glob("*.fa"))}


def screen_json(datadir: Path, sample: Path) -> dict[str, Any]:
    result = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            DB,
            "--datadir",
            str(datadir),
            "--format",
            "json",
            "--quiet",
        ],
    )
    assert result.exit_code == 0, result.stderr
    return json.loads(result.stdout)


class TestEngineSamples:
    def test_exact_locus_is_called_with_all_genes_present(
        self, cluster_db: Path, tmp_path: Path
    ) -> None:
        """Given an exact locusA contig, When screened, Then gapit.cluster/1
        names locusA best at 100/100 with all three genes present and the
        reserved phenotype null."""
        document = screen_json(cluster_db, _write_samples(cluster_db, tmp_path)["exact"])
        assert document["schema"] == "gapit.cluster/1"
        assert document["params"] == {
            "db": DB,
            "preset": "asm20",
            "min_gene_cov": 90.0,
            "min_gene_id": 90.0,
            "min_cluster_cov": 96.0,
            "threads": 1,
        }
        best = document["files"][0]["best"]
        assert best["locus"] == "locusA"
        assert best["label"] == "K antigen locus A"
        assert best["type"] == "KL101"
        assert best["coverage_pct"] == 100.0
        assert best["identity_pct"] == 100.0
        assert best["genes_present"] == 3
        assert best["genes_partial"] == 0
        assert best["genes_absent"] == 0
        assert best["phenotype"] is None
        (locus,) = document["files"][0]["loci"]
        assert locus["rank"] == 1
        assert locus["missing"] == []
        assert [gene["verdict"] for gene in locus["genes"]] == ["present"] * 3

    def test_locus_split_across_contigs_unions_to_best_call(
        self, cluster_db: Path, tmp_path: Path
    ) -> None:
        """Given locusA split 60/40 across two contigs, When screened, Then
        the union still calls locusA best with every gene present
        (fragmentation resilience)."""
        document = screen_json(cluster_db, _write_samples(cluster_db, tmp_path)["fragmented"])
        best = document["files"][0]["best"]
        assert best["locus"] == "locusA"
        assert best["coverage_pct"] == 100.0
        assert best["genes_present"] == 3
        (locus,) = document["files"][0]["loci"]
        assert {gene["gene_id"]: gene["coverage_pct"] for gene in locus["genes"]} == {
            "wzx": 100.0,
            "gtrA": 100.0,
            "manC": 100.0,
        }

    def test_truncated_locus_reports_partial_and_absent_genes(
        self, cluster_db: Path, tmp_path: Path
    ) -> None:
        """Given a contig carrying only the first 800 of 1200 locus bases,
        When screened, Then wzx is present, gtrA partial (75%), manC absent,
        and no best call is made (66.67% < the 96% cluster floor)."""
        document = screen_json(cluster_db, _write_samples(cluster_db, tmp_path)["partial"])
        entry = document["files"][0]
        assert entry["best"] is None
        (locus,) = entry["loci"]
        assert locus["locus"] == "locusA"
        assert locus["coverage_pct"] == 66.67
        verdicts = {gene["gene_id"]: gene["verdict"] for gene in locus["genes"]}
        assert verdicts == {"wzx": "present", "gtrA": "partial", "manC": "absent"}
        partial = next(gene for gene in locus["genes"] if gene["gene_id"] == "gtrA")
        assert partial["coverage_pct"] == 75.0
        assert locus["missing"] == ["gtrA", "manC"]

    def test_mutated_sample_scores_identity_near_95(self, cluster_db: Path, tmp_path: Path) -> None:
        """Given an exact locusA contig mutated at every 20th base, When
        screened, Then locus identity lands near 95% with every gene still
        present at the default identity floor."""
        document = screen_json(cluster_db, _write_samples(cluster_db, tmp_path)["mutated"])
        entry = document["files"][0]
        best = entry["best"]
        assert best is not None and best["locus"] == "locusA"
        assert 94.0 <= best["identity_pct"] <= 96.0
        (locus,) = entry["loci"]
        assert all(gene["verdict"] == "present" for gene in locus["genes"])
        assert all(94.0 <= gene["identity_pct"] <= 96.0 for gene in locus["genes"])

    def test_unrelated_sample_makes_no_call_and_no_crash(
        self, cluster_db: Path, tmp_path: Path
    ) -> None:
        """Given a contig unrelated to every locus, When screened, Then the
        loci list is empty, best is null, and the TSV renders the untyped
        row without crashing."""
        sample = _write_samples(cluster_db, tmp_path)["unrelated"]
        document = screen_json(cluster_db, sample)
        entry = document["files"][0]
        assert entry["best"] is None
        assert entry["loci"] == []
        result = runner.invoke(
            app, ["screen", str(sample), "--db", DB, "--datadir", str(cluster_db), "--quiet"]
        )
        assert result.exit_code == 0, result.stderr
        assert result.stdout.splitlines() == [
            "FILE\tBEST_LOCUS\tTYPE\tCOVERAGE\tIDENTITY\tPRESENT\tPARTIAL\tMISSING_IDS",
            f"{sample}\t-\t-\t0.00\t0.00\t0\t0\t-",
        ]

    def test_min_cluster_cov_floor_gates_the_best_call(
        self, cluster_db: Path, tmp_path: Path
    ) -> None:
        """Given the mutated sample (99.92% locus coverage) and
        --min-cluster-cov 99.95, When screened, Then no best call is made
        even though locusA ranks first with covered loci."""
        sample = _write_samples(cluster_db, tmp_path)["mutated"]
        result = runner.invoke(
            app,
            [
                "screen",
                str(sample),
                "--db",
                DB,
                "--datadir",
                str(cluster_db),
                "--format",
                "json",
                "--quiet",
                "--min-cluster-cov",
                "99.95",
            ],
        )
        assert result.exit_code == 0, result.stderr
        entry = json.loads(result.stdout)["files"][0]
        assert entry["best"] is None
        assert entry["loci"] and entry["loci"][0]["locus"] == "locusA"

    def test_debug_echoes_asm20_argv(self, cluster_db: Path, tmp_path: Path) -> None:
        """Given --debug, When screening, Then the minimap2 argv line on
        stderr carries -x asm20 and --cs (engine telemetry), and stdout is
        unchanged from the plain run."""
        sample = _write_samples(cluster_db, tmp_path)["exact"]
        args = ["screen", str(sample), "--db", DB, "--datadir", str(cluster_db), "--format", "json"]
        debug_run = runner.invoke(app, [*args, "--debug", "--quiet"])
        plain_run = runner.invoke(app, [*args, "--quiet"])
        assert debug_run.exit_code == 0, debug_run.stderr
        argv_lines = [
            line for line in debug_run.stderr.splitlines() if line.startswith("gapit: run:")
        ]
        assert argv_lines and "-x asm20" in argv_lines[0] and "--cs" in argv_lines[0]
        assert debug_run.stdout == plain_run.stdout


class TestFormats:
    def test_tsv_default_header_and_multiple_files(self, cluster_db: Path, tmp_path: Path) -> None:
        """Given two samples, When screened with the default TSV output,
        Then one row per file follows the cluster header in input order."""
        samples = _write_samples(cluster_db, tmp_path)
        result = runner.invoke(
            app,
            [
                "screen",
                str(samples["exact"]),
                str(samples["unrelated"]),
                "--db",
                DB,
                "--datadir",
                str(cluster_db),
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.stderr
        header, exact_row, unrelated_row = result.stdout.splitlines()
        assert header == "FILE\tBEST_LOCUS\tTYPE\tCOVERAGE\tIDENTITY\tPRESENT\tPARTIAL\tMISSING_IDS"
        assert exact_row.startswith(f"{samples['exact']}\tlocusA\tKL101\t100.00\t100.00\t3\t0\t-")
        assert unrelated_row == f"{samples['unrelated']}\t-\t-\t0.00\t0.00\t0\t0\t-"

    def test_csv_and_noheader_nopath(self, cluster_db: Path, tmp_path: Path) -> None:
        """Given --format csv --noheader --nopath, When screened, Then the
        header is gone, fields are comma-joined, and FILE is a basename."""
        sample = _write_samples(cluster_db, tmp_path)["exact"]
        result = runner.invoke(
            app,
            [
                "screen",
                str(sample),
                "--db",
                DB,
                "--datadir",
                str(cluster_db),
                "--quiet",
                "--format",
                "csv",
                "--noheader",
                "--nopath",
            ],
        )
        assert result.exit_code == 0, result.stderr
        (row,) = result.stdout.splitlines()
        assert row.startswith("exact.fa,locusA,KL101,100.00,100.00,3,0,-")

    def test_md_frontmatter_and_gene_table(self, cluster_db: Path, tmp_path: Path) -> None:
        """Given --format md, When screened, Then YAML frontmatter carries
        the engine params and the best locus's gene table shows verdicts."""
        sample = _write_samples(cluster_db, tmp_path)["exact"]
        result = runner.invoke(
            app,
            [
                "screen",
                str(sample),
                "--db",
                DB,
                "--datadir",
                str(cluster_db),
                "--format",
                "md",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.stderr
        text = result.stdout
        assert text.startswith("---\n")
        assert "schema: gapit.cluster/1" in text
        assert "preset: asm20" in text
        assert "min_cluster_cov: 96.0" in text
        assert "## `ctg_exact`" not in text  # sections key on the FILE, not the contig
        assert f"## `{sample}`" in text
        for column in ("Gene", "Verdict", "wzx", "manC", "present"):
            assert column in text


class TestGoldens:
    """Byte-exact renders at a pinned timestamp (minimap2 2.31-r1302)."""

    @staticmethod
    def _report(datadir: Path, sample: Path) -> Any:
        (database,) = discover_databases(datadir)
        return screen_cluster_file(sample, database, load_features(database), PARAMS)

    def _check(self, name: str, rendered: str) -> None:
        golden = GOLDEN / f"cluster_{name}"
        assert rendered == golden.read_text(encoding="utf-8"), f"golden mismatch: {golden}"

    @pytest.mark.parametrize(
        ("name", "sample_name"), [("exact", "exact"), ("fragmented", "fragmented")]
    )
    def test_cluster_goldens(
        self,
        cluster_db: Path,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        name: str,
        sample_name: str,
    ) -> None:
        """Given the exact/fragmented samples, When the engine renders all
        three formats at the pinned now, Then each is byte-identical to the
        committed golden."""
        sample = _write_samples(cluster_db, tmp_path)[sample_name]
        monkeypatch.chdir(sample.parent)
        relative = Path(sample.name)
        report = self._report(cluster_db, relative)
        self._check(f"{name}.json", render_cluster_json([report], PARAMS, now=PINNED_NOW))
        self._check(f"{name}.md", render_cluster_md([report], PARAMS, now=PINNED_NOW))
        self._check(
            f"{name}.tsv", format_cluster_tsv([report], csv=False, noheader=False, nopath=False)
        )
        parsed = json.loads(render_cluster_json([report], PARAMS, now=PINNED_NOW))
        assert parsed["schema"] == "gapit.cluster/1"


class TestKindGuards:
    """The shared db-kind guard, both directions, plus reads-mode."""

    def _usage(self, result: Result, message_part: str) -> None:
        assert result.exit_code == 2
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "USAGE_ERROR"
        assert message_part in envelope["message"]

    @pytest.mark.parametrize(
        ("extra", "message_part"),
        [
            (["--minid", "85"], "--min-gene-cov/--min-gene-id"),
            (["--mincov", "85"], "--min-gene-cov/--min-gene-id"),
            (["--merge-fragments"], "--merge-fragments is not available with a cluster database"),
            (["--jobs", "2"], "--jobs is not available with a cluster database"),
            (["--aligner", "blastn"], "--aligner is not available with a cluster database"),
            (["--min-gene-cov", "150"], "--min-gene-cov must be in [0, 100]"),
        ],
    )
    def test_gene_flags_on_cluster_db_are_usage_errors(
        self, cluster_db: Path, tmp_path: Path, extra: list[str], message_part: str
    ) -> None:
        """Given a cluster db and any gene-engine flag at a non-default
        value, When screened, Then exit 2 with a message naming the cluster
        flags."""
        sample = _write_samples(cluster_db, tmp_path)["exact"]
        result = runner.invoke(
            app,
            ["screen", str(sample), "--db", DB, "--datadir", str(cluster_db), "--quiet", *extra],
        )
        self._usage(result, message_part)

    def test_cluster_flags_on_gene_db_are_usage_error(self, tmp_path: Path) -> None:
        """Given a gene db and a cluster threshold flag, When screened,
        Then exit 2 naming the cluster-only flags."""
        datadir = tmp_path / "datadir"
        shutil.copytree(GENE_DB, datadir)
        contig = tmp_path / "contigs.fna"
        shutil.copyfile(Path(__file__).parent / "data" / "contigs" / "full.fa", contig)
        for flag in ("--min-gene-cov", "--min-gene-id", "--min-cluster-cov"):
            value = "95" if flag != "--min-cluster-cov" else "97"
            result = runner.invoke(
                app,
                ["screen", str(contig), "--db", "tinyamr", "--datadir", str(datadir), flag, value],
            )
            self._usage(result, "--min-gene-cov/--min-gene-id/--min-cluster-cov apply to cluster")

    def test_reads_mode_rejects_cluster_db(self, cluster_db: Path, tmp_path: Path) -> None:
        """Given --r1 reads with a cluster db, When screened, Then exit 2:
        cluster screening is assembly-only in v1."""
        reads = Path(__file__).parent / "data" / "reads" / "tetx_full.fq"
        result = runner.invoke(
            app, ["screen", "--r1", str(reads), "--db", DB, "--datadir", str(cluster_db)]
        )
        self._usage(result, "cluster databases are assembly-contig screening only")

    def test_aligner_minimap2_rejects_cluster_db(self, cluster_db: Path, tmp_path: Path) -> None:
        """Given --aligner minimap2 with a cluster db, When screened, Then
        exit 2 (the reads pipeline never screens cluster dbs)."""
        sample = _write_samples(cluster_db, tmp_path)["exact"]
        result = runner.invoke(
            app,
            [
                "screen",
                str(sample),
                "--db",
                DB,
                "--datadir",
                str(cluster_db),
                "--aligner",
                "minimap2",
            ],
        )
        self._usage(result, "cluster databases are assembly-contig screening only")

    def test_reads_mode_rejects_cluster_flags(self, cluster_db: Path, tmp_path: Path) -> None:
        """Given --r1 with a cluster threshold flag, When screened, Then
        exit 2 — cluster flags never ride along in reads mode."""
        reads = Path(__file__).parent / "data" / "reads" / "tetx_full.fq"
        gene_dir = tmp_path / "gene"
        shutil.copytree(GENE_DB, gene_dir)
        result = runner.invoke(
            app,
            [
                "screen",
                "--r1",
                str(reads),
                "--db",
                "tinyamr",
                "--datadir",
                str(gene_dir),
                "--min-gene-cov",
                "95",
            ],
        )
        self._usage(result, "--min-gene-cov/--min-gene-id/--min-cluster-cov apply to cluster")


class TestSchema:
    def test_gapit_schema_cluster(self) -> None:
        """Given the schema command, When asked for cluster, Then the
        gapit.cluster/1 JSON Schema prints with the reserved phenotype."""
        result = runner.invoke(app, ["schema", "cluster"])
        assert result.exit_code == 0
        payload = json.loads(result.stdout)
        assert payload["title"] == "ClusterDocument"
        assert "phenotype" in payload["$defs"]["ClusterBestDocument"]["properties"]


class TestMcp:
    def _call(self, name: str, arguments: dict[str, object]) -> tuple[bool, str]:
        raw = (
            json.dumps(
                {
                    "jsonrpc": "2.0",
                    "id": 1,
                    "method": "tools/call",
                    "params": {"name": name, "arguments": arguments},
                }
            )
            + "\n"
        )
        out = StringIO()
        serve(StringIO(raw), out)
        (response,) = (json.loads(line) for line in out.getvalue().splitlines())
        result = response["result"]
        return result["isError"], result["content"][0]["text"]

    def test_mcp_screen_cluster_db_returns_cluster_document(
        self, cluster_db: Path, tmp_path: Path
    ) -> None:
        """Given the MCP screen tool and a cluster-kind db, When called with
        assembly files, Then it returns the gapit.cluster/1 document with the
        best-locus call; minGeneId overrides flow through."""
        sample = _write_samples(cluster_db, tmp_path)["mutated"]
        is_error, text = self._call(
            "screen", {"files": [str(sample)], "db": DB, "datadir": str(cluster_db)}
        )
        assert not is_error
        document = json.loads(text)
        assert document["schema"] == "gapit.cluster/1"
        assert document["files"][0]["best"]["locus"] == "locusA"

        stricter, stricter_text = self._call(
            "screen",
            {"files": [str(sample)], "db": DB, "datadir": str(cluster_db), "minGeneId": 96},
        )
        assert not stricter
        entry = json.loads(stricter_text)["files"][0]
        assert entry["best"]["genes_partial"] == 3
        assert entry["loci"][0]["missing"] == ["wzx", "gtrA", "manC"]

    def test_mcp_gene_params_rejected_on_cluster_db(self, cluster_db: Path, tmp_path: Path) -> None:
        """Given the MCP screen tool, a cluster db, and the blastn-only
        minid, When called, Then isError carries the gapit.error/1 usage
        envelope naming the cluster flags."""
        sample = _write_samples(cluster_db, tmp_path)["exact"]
        is_error, text = self._call(
            "screen", {"files": [str(sample)], "db": DB, "datadir": str(cluster_db), "minid": 90}
        )
        assert is_error
        assert json.loads(text)["code"] == "USAGE_ERROR"


@pytest.mark.skipif(not REAL_CPS.is_file(), reason=f"real CPS locus gbk not present: {REAL_CPS}")
class TestRealDataSmoke:
    def test_cpsgc_locus_ranks_first(self, tmp_path: Path) -> None:
        """Given a db built from the real V. parahaemolyticus CPSgc.gbk and
        a synthetic contig extracted from one of its loci, When screened,
        Then that locus ranks first (not committed; smoke only)."""
        datadir = tmp_path / "datadir"
        datadir.mkdir()
        build = runner.invoke(
            app, ["db", "build", "cpsreal", str(REAL_CPS), "--datadir", str(datadir)]
        )
        assert build.exit_code == 0, build.stderr
        sequences = {
            record.id: record.sequence for record in iter_fasta(datadir / "cpsreal" / "sequences")
        }
        locus_id, locus_sequence = next(iter(sequences.items()))
        contig = tmp_path / "smoke.fa"
        contig.write_text(f">smoke_{locus_id}\n{locus_sequence}\n", encoding="utf-8")

        result = runner.invoke(
            app,
            [
                "screen",
                str(contig),
                "--db",
                "cpsreal",
                "--datadir",
                str(datadir),
                "--format",
                "json",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.stderr
        entry = json.loads(result.stdout)["files"][0]
        assert entry["loci"], "expected at least one covered locus"
        assert entry["loci"][0]["rank"] == 1
        assert entry["loci"][0]["locus"] == locus_id
        assert entry["loci"][0]["coverage_pct"] == 100.0
        assert entry["best"]["locus"] == locus_id
