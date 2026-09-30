"""Cluster-typing integration suites (stage 3): a typed fixture database
(two weighted_genes rules / two cluster_match rules over tests/data/cluster/
screening.gbk) screened end to end through the real engine.

Covers the output wiring (phenotype + phenotype_detail in gapit.cluster/1,
the PHENOTYPE TSV column variant, the Markdown phenotype line), the
ambiguity decision on a two-locus sample, goldens for both rule kinds at a
pinned timestamp, and the eval-time TYPING_UNKNOWN_GENE guard. The stage-2
untyped goldens stay byte-identical (their own suite); the pure scoring
math lives in test_typing_engine.py.
"""

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from gapit.cli import app
from gapit.cluster import ClusterParams, load_features, load_typing, screen_cluster_file
from gapit.db import discover_databases
from gapit.fasta import iter_fasta
from gapit.formats.cluster import format_cluster_tsv, render_cluster_json
from gapit.formats.cluster_md import render_cluster_md
from gapit.typing_engine import evaluate_typing

DATA = Path(__file__).parent / "data" / "cluster"
GOLDEN = Path(__file__).parent / "golden"
DB = "cps"
PINNED_NOW = datetime(2026, 9, 30, 12, 0, 0, tzinfo=UTC)
PARAMS = ClusterParams(db=DB)

runner = CliRunner()


@pytest.fixture()
def typed_db(tmp_path: Path) -> Path:
    """Build the screening fixture db with the weighted_genes typing doc
    into a fresh datadir (the cluster_match variant is built per test via
    _build_typed)."""
    return _build_typed(tmp_path, "typing_screen.json", tag="fixture")


@pytest.fixture()
def samples(tmp_path: Path, typed_db: Path) -> dict[str, Path]:
    """exact = locusA alone; both = locusA + locusB on two contigs."""
    loci = {record.id: record.sequence for record in iter_fasta(typed_db / DB / "sequences")}
    out = tmp_path / "samples"
    out.mkdir()
    (out / "exact.fa").write_text(f">ctg_a\n{loci['locusA']}\n", encoding="utf-8")
    (out / "both.fa").write_text(
        f">ctg_a\n{loci['locusA']}\n>ctg_b\n{loci['locusB']}\n", encoding="utf-8"
    )
    return {path.stem: path for path in sorted(out.glob("*.fa"))}


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


class TestOutputWiring:
    def test_weighted_genes_phenotype_with_breakdown(
        self, typed_db: Path, samples: dict[str, Path]
    ) -> None:
        """Given the exact locusA sample against the weighted_genes db,
        When screened as JSON, Then best carries phenotype K101 with a
        detail whose components are the per-gene contributions and whose
        runner-up is K102."""
        document = screen_json(typed_db, samples["exact"])
        best = document["files"][0]["best"]
        assert best["locus"] == "locusA"
        assert best["phenotype"] == "K101"
        detail = best["phenotype_detail"]
        assert detail["score"] == 1.0
        assert detail["confidence"] == "high"
        assert detail["components"] == [
            {"name": "wzx", "score": 0.3333},
            {"name": "gtrA", "score": 0.3333},
            {"name": "manC", "score": 0.3333},
        ]
        assert detail["runner_up"] == {"phenotype": "K102", "score": 0.0}
        assert detail["ambiguous"] == []

    def test_cluster_match_phenotype_with_components(
        self, tmp_path: Path, samples: dict[str, Path]
    ) -> None:
        """Given the exact locusA sample against the cluster_match db, When
        screened, Then the detail components are coverage/identity/key_genes
        with their weighted contributions."""
        datadir = _build_typed(tmp_path, "typing_screen_cm.json")
        document = screen_json(datadir, samples["exact"])
        detail = document["files"][0]["best"]["phenotype_detail"]
        assert detail["score"] == 1.0
        assert [component["name"] for component in detail["components"]] == [
            "coverage",
            "identity",
            "key_genes",
        ]
        assert [component["score"] for component in detail["components"]] == [0.5, 0.3, 0.2]

    def test_tsv_gains_phenotype_column_only_when_typed(
        self, typed_db: Path, samples: dict[str, Path], tmp_path: Path
    ) -> None:
        """Given the typed db, When screened as TSV, Then the header gains
        PHENOTYPE after TYPE and the row carries the called phenotype; an
        untyped rebuild of the same fixture keeps the stage-2 header."""
        result = runner.invoke(
            app,
            ["screen", str(samples["exact"]), "--db", DB, "--datadir", str(typed_db), "--quiet"],
        )
        assert result.exit_code == 0, result.stderr
        header, row = result.stdout.splitlines()
        assert header == (
            "FILE\tBEST_LOCUS\tTYPE\tPHENOTYPE\tCOVERAGE\tIDENTITY\tPRESENT\tPARTIAL\tMISSING_IDS"
        )
        assert "\tlocusA\tKL101\t" in row
        assert "\tK101\t100.00\t100.00\t" in row

        untyped = tmp_path / "untyped"
        untyped.mkdir()
        build = runner.invoke(
            app, ["db", "build", DB, str(DATA / "screening.gbk"), "--datadir", str(untyped)]
        )
        assert build.exit_code == 0, build.stderr
        plain = runner.invoke(
            app,
            ["screen", str(samples["exact"]), "--db", DB, "--datadir", str(untyped), "--quiet"],
        )
        assert plain.exit_code == 0, plain.stderr
        assert plain.stdout.splitlines()[0] == (
            "FILE\tBEST_LOCUS\tTYPE\tCOVERAGE\tIDENTITY\tPRESENT\tPARTIAL\tMISSING_IDS"
        )

    def test_md_carries_phenotype_line_and_column(
        self, tmp_path: Path, samples: dict[str, Path]
    ) -> None:
        """Given the typed db, When screened as Markdown, Then the summary
        table carries the Phenotype column and the file section opens with
        the phenotype score line."""
        datadir = _build_typed(tmp_path, "typing_screen.json")
        result = runner.invoke(
            app,
            [
                "screen",
                str(samples["exact"]),
                "--db",
                DB,
                "--datadir",
                str(datadir),
                "--format",
                "md",
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.stderr
        assert "| Phenotype |" in result.stdout
        assert "| locusA | KL101 | K101 |" in result.stdout
        assert "Phenotype `K101` (score 1.0000, high confidence)" in result.stdout

    def test_mcp_screen_returns_typed_cluster_document(
        self, typed_db: Path, samples: dict[str, Path]
    ) -> None:
        """Given the MCP screen tool and the typed db, When called, Then
        the cluster document carries the phenotype call (shared use-case,
        same JSON as the CLI)."""
        from io import StringIO

        from gapit.mcp import serve

        raw = (
            json.dumps(
                {
                    "jsonrpc": "2.0",
                    "id": 1,
                    "method": "tools/call",
                    "params": {
                        "name": "screen",
                        "arguments": {
                            "files": [str(samples["exact"])],
                            "db": DB,
                            "datadir": str(typed_db),
                        },
                    },
                }
            )
            + "\n"
        )
        out = StringIO()
        serve(StringIO(raw), out)
        (response,) = (json.loads(line) for line in out.getvalue().splitlines())
        assert not response["result"]["isError"]
        best = json.loads(response["result"]["content"][0]["text"])["files"][0]["best"]
        assert best["phenotype"] == "K101"
        assert best["phenotype_detail"]["confidence"] == "high"


class TestAmbiguity:
    def test_two_matching_loci_yield_ambiguous_call(
        self, tmp_path: Path, samples: dict[str, Path]
    ) -> None:
        """Given a sample covering locusA AND locusB fully, When screened
        against the weighted_genes db (K101 and K102 both fire), Then the
        phenotype is null with confidence ambiguous and the two candidates
        are listed in order; the TSV phenotype cell is '-'."""
        datadir = _build_typed(tmp_path, "typing_screen.json")
        document = screen_json(datadir, samples["both"])
        best = document["files"][0]["best"]
        assert best["locus"] == "locusA"
        assert best["phenotype"] is None
        detail = best["phenotype_detail"]
        assert detail["confidence"] == "ambiguous"
        assert detail["ambiguous"] == [
            {"phenotype": "K101", "score": 1.0},
            {"phenotype": "K102", "score": 1.0},
        ]
        result = runner.invoke(
            app,
            ["screen", str(samples["both"]), "--db", DB, "--datadir", str(datadir), "--quiet"],
        )
        assert result.exit_code == 0, result.stderr
        row = result.stdout.splitlines()[1]
        assert row.split("\t")[3] == "-"

    def test_unrelated_sample_renders_untyped_dash_row(
        self, tmp_path: Path, typed_db: Path
    ) -> None:
        """Given a sample with no locus coverage, When screened against the
        typed db, Then the TSV row is the no-call row with a '-' phenotype
        and no phenotype_detail appears anywhere."""
        sample = tmp_path / "unrelated.fa"
        sample.write_text(">ctg\n" + "ACGT" * 250 + "\n", encoding="utf-8")
        document = screen_json(typed_db, sample)
        assert document["files"][0]["best"] is None
        result = runner.invoke(
            app, ["screen", str(sample), "--db", DB, "--datadir", str(typed_db), "--quiet"]
        )
        assert result.exit_code == 0, result.stderr
        row = result.stdout.splitlines()[1]
        assert row.endswith("\t-\t-\t-\t0.00\t0.00\t0\t0\t-")


class TestGoldens:
    """Byte-exact typed renders at a pinned timestamp, both rule kinds."""

    @staticmethod
    def _report(datadir: Path, sample: Path) -> Any:
        (database,) = discover_databases(datadir)
        features = load_features(database)
        typing_document = load_typing(database)
        assert typing_document is not None
        report = screen_cluster_file(sample, database, features, PARAMS)
        return evaluate_typing(report, typing_document)

    def _check(self, name: str, rendered: str) -> None:
        golden = GOLDEN / f"cluster_{name}"
        assert rendered == golden.read_text(encoding="utf-8"), f"golden mismatch: {golden}"

    @pytest.mark.parametrize(
        ("fixture", "sample", "prefix"),
        [
            ("typing_screen.json", "exact", "typing_exact"),
            ("typing_screen_cm.json", "exact", "typing_cm_exact"),
            ("typing_screen.json", "both", "typing_ambiguous"),
        ],
    )
    def test_typed_goldens(
        self,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        fixture: str,
        sample: str,
        prefix: str,
    ) -> None:
        """Given the typed fixture dbs, When the engine renders all three
        formats at the pinned now, Then each is byte-identical to the
        committed golden."""
        datadir = _build_typed(tmp_path, fixture)
        loci = {record.id: record.sequence for record in iter_fasta(datadir / DB / "sequences")}
        out = tmp_path / "samples"
        out.mkdir()
        (out / "exact.fa").write_text(f">ctg_a\n{loci['locusA']}\n", encoding="utf-8")
        (out / "both.fa").write_text(
            f">ctg_a\n{loci['locusA']}\n>ctg_b\n{loci['locusB']}\n", encoding="utf-8"
        )
        sample_path = out / f"{sample}.fa"
        monkeypatch.chdir(out)
        relative = Path(sample_path.name)
        report = self._report(datadir, relative)
        self._check(
            f"{prefix}.json", render_cluster_json([report], PARAMS, now=PINNED_NOW, typed=True)
        )
        self._check(f"{prefix}.md", render_cluster_md([report], PARAMS, now=PINNED_NOW, typed=True))
        self._check(
            f"{prefix}.tsv",
            format_cluster_tsv([report], csv=False, noheader=False, nopath=False, typed=True),
        )


class TestEvalValidation:
    def test_unknown_gene_aborts_screen_with_typed_error(
        self, tmp_path: Path, samples: dict[str, Path]
    ) -> None:
        """Given a db whose typing.json references a gene no locus carries
        (installed after the build), When screened, Then exit 4 with
        TYPING_UNKNOWN_GENE naming the gene — before any minimap2 run."""
        datadir = _build_typed(tmp_path, "typing_screen.json")
        (datadir / DB / "typing.json").write_text(
            json.dumps(
                {
                    "schema": "gapit.typing/1",
                    "rules": [
                        {
                            "model": "weighted_genes",
                            "phenotype": "KX",
                            "weights": {"ghost_gene": 1.0},
                            "identity_floor": 90.0,
                        }
                    ],
                    "cutoff": 0.9,
                    "ambiguity_margin": 0.05,
                    "fallback": "unknown",
                }
            ),
            encoding="utf-8",
        )
        result = runner.invoke(
            app, ["screen", str(samples["exact"]), "--db", DB, "--datadir", str(datadir), "--quiet"]
        )
        assert result.exit_code == 4
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "TYPING_UNKNOWN_GENE"
        assert envelope["context"]["gene"] == "ghost_gene"
        assert "gapit: run:" not in result.stderr  # failed before any minimap2 invocation


def _build_typed(tmp_path: Path, fixture: str, *, tag: str | None = None) -> Path:
    """Build the screening fixture db with one typing spec (test helper;
    one datadir per fixture/tag so combined-fixture tests never collide)."""
    datadir = tmp_path / f"datadir-{tag or fixture}"
    datadir.mkdir()
    result = runner.invoke(
        app,
        [
            "db",
            "build",
            DB,
            str(DATA / "screening.gbk"),
            "--datadir",
            str(datadir),
            "--typing",
            str(DATA / fixture),
        ],
    )
    assert result.exit_code == 0, result.stderr
    return datadir
