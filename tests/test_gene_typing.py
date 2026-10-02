"""Gene-path typing integration suites (typing/2 stage 1): a typed marker
fixture database (four marker genes, two weighted_genes schemes — DEC
pathology-panel style) built with ``db build --typing`` and run through
the two-stage pipeline: ``screen -o result.tsv`` (real blastn engine),
then ``gapit typing result.tsv``.

Covers the gapit.typing_result/1 calls (JSON + MD + TSV), the
identity-floor gating over the folded best rows, the multi-scheme
independence (one scheme ambiguous while another calls), fallback calls
on no-hit files, the v1-document degradation to the ``default`` scheme
key, the typed/untyped screen byte-identity (the purification pin: gene
screening is pure detection), the reads-mode availability of typed gene
dbs, the build-time TYPING_UNKNOWN_GENE and rule-kind validation, and
goldens at a pinned timestamp. The cluster-side typing suites live in
test_cluster_typing.py; the pure scoring math in test_typing_engine.py;
the command's input parsing in test_typing_command.py.
"""

import json
import re
from datetime import UTC, datetime
from io import StringIO
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from gapit.blast import ensure_blast, screen_file
from gapit.cli import app
from gapit.db import discover_databases
from gapit.formats.json import render_json
from gapit.formats.md import render_markdown
from gapit.mcp import serve
from gapit.report import ScreeningParams

DATA = Path(__file__).parent / "data" / "typing"
GOLDEN = Path(__file__).parent / "golden"
DB = "markers"
PINNED_NOW = datetime(2026, 9, 30, 12, 0, 0, tzinfo=UTC)
PARAMS = ScreeningParams(db=DB)

runner = CliRunner()

SAMPLES = ("exact", "partial", "mutated", "negative")


@pytest.fixture(scope="module")
def typed_db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The marker fixture db with the two-scheme typing/2 document, plus an
    untyped rebuild of the same FASTA in a sibling datadir."""
    datadir = tmp_path_factory.mktemp("typed")
    result = runner.invoke(
        app,
        [
            "db",
            "build",
            DB,
            str(DATA / "markers.fa"),
            "--datadir",
            str(datadir),
            "--typing",
            str(DATA / "typing_markers_v2.json"),
        ],
    )
    assert result.exit_code == 0, result.stderr
    return datadir


@pytest.fixture(scope="module")
def untyped_db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    datadir = tmp_path_factory.mktemp("untyped")
    result = runner.invoke(
        app, ["db", "build", DB, str(DATA / "markers.fa"), "--datadir", str(datadir)]
    )
    assert result.exit_code == 0, result.stderr
    return datadir


def screen_json(datadir: Path, sample: str) -> dict[str, Any]:
    result = runner.invoke(
        app,
        [
            "screen",
            str(DATA / f"{sample}.fa"),
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


def screen_table(datadir: Path, samples: list[str], tmp_path: Path) -> Path:
    """Stage 1 of the pipeline: screen the sample(s) into one TSV table."""
    table = tmp_path / "screen.tsv"
    result = runner.invoke(
        app,
        [
            "screen",
            *(str(DATA / f"{sample}.fa") for sample in samples),
            "--db",
            DB,
            "--datadir",
            str(datadir),
            "-o",
            str(table),
            "--quiet",
        ],
    )
    assert result.exit_code == 0, result.stderr
    return table


def typing_json(datadir: Path, sample: str, tmp_path: Path) -> dict[str, Any]:
    """Stage 2: type the screened table, returning the typing_result/1
    document."""
    result = runner.invoke(
        app,
        [
            "typing",
            str(screen_table(datadir, [sample], tmp_path)),
            "--datadir",
            str(datadir),
            "--format",
            "json",
            "--quiet",
        ],
    )
    assert result.exit_code == 0, result.stderr
    return json.loads(result.stdout)


def typing_stdout(datadir: Path, samples: list[str], tmp_path: Path, *extra: str) -> str:
    result = runner.invoke(
        app,
        [
            "typing",
            str(screen_table(datadir, samples, tmp_path)),
            "--datadir",
            str(datadir),
            "--quiet",
            *extra,
        ],
    )
    assert result.exit_code == 0, result.stderr
    return result.stdout


class TestPhenotypeCalls:
    def test_exact_sample_calls_each_scheme(self, typed_db: Path, tmp_path: Path) -> None:
        """Given the exact sample (every marker at 100%), When typed, Then
        pathotype calls EHEC (the marker_c negative weight suppresses EPEC)
        while toxin is ambiguous on the Toxin1/Toxin2 tie — the schemes
        decide independently."""
        phenotypes = typing_json(typed_db, "exact", tmp_path)["files"][0]["phenotypes"]
        assert phenotypes["pathotype"]["phenotype"] == "EHEC"
        assert phenotypes["pathotype"]["confidence"] == "high"
        assert phenotypes["pathotype"]["score"] == 1.0
        assert phenotypes["pathotype"]["runner_up"] == {"phenotype": "EPEC", "score": 0.0}
        assert phenotypes["pathotype"]["components"] == [
            {"name": "marker_a", "score": 0.5},
            {"name": "marker_c", "score": 0.5},
        ]
        assert phenotypes["toxin"]["phenotype"] is None
        assert phenotypes["toxin"]["confidence"] == "ambiguous"
        assert phenotypes["toxin"]["ambiguous"] == [
            {"phenotype": "Toxin1", "score": 1.0},
            {"phenotype": "Toxin2", "score": 1.0},
        ]

    def test_partial_sample_calls_the_matching_rules(self, typed_db: Path, tmp_path: Path) -> None:
        """Given the partial sample (marker_a + marker_b only), When typed,
        Then pathotype calls EPEC (EHEC's require_any on marker_c is unmet,
        and with marker_c absent nothing suppresses EPEC) and toxin calls
        Toxin1."""
        phenotypes = typing_json(typed_db, "partial", tmp_path)["files"][0]["phenotypes"]
        assert phenotypes["pathotype"]["phenotype"] == "EPEC"
        assert phenotypes["pathotype"]["confidence"] == "high"
        assert phenotypes["toxin"]["phenotype"] == "Toxin1"

    def test_identity_floor_gates_a_mutated_marker(self, typed_db: Path, tmp_path: Path) -> None:
        """Given the mutated sample (marker_b at ~93.9% identity — above the
        run's minid but below Toxin1's 95 identity_floor), When typed, Then
        the hit survives into the table yet Toxin1 scores 0 and Toxin2 is
        called: floor semantics are the cluster path's."""
        document = screen_json(typed_db, "mutated")
        hits = document["files"][0]["hits"]
        marker_b = next(hit for hit in hits if hit["gene"] == "marker_b")
        assert 93.0 < marker_b["identity_pct"] < 95.0
        phenotypes = typing_json(typed_db, "mutated", tmp_path)["files"][0]["phenotypes"]
        assert phenotypes["toxin"]["phenotype"] == "Toxin2"
        assert phenotypes["toxin"]["components"] == [{"name": "marker_d", "score": 1.0}]
        assert phenotypes["toxin"]["runner_up"] == {"phenotype": "Toxin1", "score": 0.0}
        assert phenotypes["pathotype"]["phenotype"] == "EHEC"

    def test_negative_sample_falls_back_per_scheme(self, typed_db: Path, tmp_path: Path) -> None:
        """Given a screened sample with no hits at all, When its table types,
        Then the file is absent from the output — presence is row
        existence, and a table of only no-hit files is the TYPING_NO_DATA
        input error (covered in test_typing_command.py)."""
        document = json.loads(typing_stdout(typed_db, list(SAMPLES), tmp_path, "--format", "json"))
        files = [entry["file"] for entry in document["files"]]
        assert str(DATA / "negative.fa") not in files
        assert str(DATA / "exact.fa") in files


class TestPurifiedScreen:
    @staticmethod
    def _without_created_at(document: str) -> str:
        return re.sub(r"created_at: [^\\\n]+|\"created_at\": \"[^\"]+\"", "created_at: -", document)

    def test_typed_and_untyped_screens_are_byte_identical(
        self, typed_db: Path, untyped_db: Path
    ) -> None:
        """Given the same inputs against the typed and the untyped build,
        When screened in every format, Then the outputs are byte-identical
        (modulo each run's created_at timestamp) — gene screening is pure
        detection; designation lives in `gapit typing` (the inverted former
        phenotype-column pin)."""
        inputs = [str(DATA / f"{sample}.fa") for sample in SAMPLES]
        for extra in ([], ["--format", "json"], ["--format", "md"]):
            typed = runner.invoke(
                app, ["screen", *inputs, "--db", DB, "--datadir", str(typed_db), "--quiet", *extra]
            )
            untyped = runner.invoke(
                app,
                ["screen", *inputs, "--db", DB, "--datadir", str(untyped_db), "--quiet", *extra],
            )
            assert typed.exit_code == 0, typed.stderr
            assert untyped.exit_code == 0, untyped.stderr
            left, right = map(self._without_created_at, (typed.stdout, untyped.stdout))
            assert left == right
            assert "phenotypes" not in typed.stdout
            assert "Phenotype" not in typed.stdout

    def test_untyped_json_has_no_phenotypes_anywhere(self, untyped_db: Path) -> None:
        """Given the untyped rebuild, When screened as JSON, Then no
        ``phenotypes`` key appears on any file block."""
        for sample in SAMPLES:
            document = screen_json(untyped_db, sample)
            for entry in document["files"]:
                assert "phenotypes" not in entry

    def test_aligner_minimap2_screens_a_typed_gene_db(self, typed_db: Path) -> None:
        """Given the typed gene db and --aligner minimap2, When screening,
        Then it succeeds — the former contig-mode-only guard existed to
        protect inline typing annotations, which no longer ride the
        screen."""
        result = runner.invoke(
            app,
            [
                "screen",
                str(DATA / "exact.fa"),
                "--aligner",
                "minimap2",
                "--db",
                DB,
                "--datadir",
                str(typed_db),
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.stderr


class TestTypingOutput:
    def test_tsv_renders_one_row_per_scheme_with_the_pair_in_notes(
        self, typed_db: Path, tmp_path: Path
    ) -> None:
        """Given the exact sample whose toxin scheme is ambiguous, When the
        table types as TSV, Then PHENOTYPE renders ``-`` for that scheme
        with the candidate pair carried in NOTES, while called schemes
        render their phenotype and runner-up; the GENES cell after PHENOTYPE repeats
        the FILE's sorted ``;``-joined gene list on every scheme row."""
        output = typing_stdout(typed_db, ["exact"], tmp_path)
        lines = output.splitlines()
        genes = "marker_a;marker_b;marker_c;marker_d"
        assert lines[0] == "FILE\tSCHEME\tPHENOTYPE\tGENES\tCONFIDENCE\tSCORE\tRUNNER_UP\tNOTES"
        pathotype, toxin = lines[1], lines[2]
        assert pathotype == (
            f"{DATA / 'exact.fa'}\tpathotype\tEHEC\t{genes}\thigh\t1.0000\tEPEC (0.0000)\t"
        )
        assert toxin == (
            f"{DATA / 'exact.fa'}\ttoxin\t-\t{genes}\tambiguous\t1.0000\t-\t"
            f"ambiguous: Toxin1 (1.0000), Toxin2 (1.0000)"
        )

    def test_md_renders_the_eight_column_table(self, typed_db: Path, tmp_path: Path) -> None:
        """Given the typed db and the exact sample, When typed as Markdown,
        Then one table row per scheme appears (the ambiguous toxin row
        carrying the pair in NOTES, the GENES cell following PHENOTYPE),
        frontmatter naming the db and the source table."""
        output = typing_stdout(typed_db, ["exact"], tmp_path, "--format", "md")
        assert "schema: gapit.typing_result/1" in output
        assert f"db: {DB}" in output
        assert "source:" in output
        header = "| FILE | SCHEME | PHENOTYPE | GENES | CONFIDENCE | SCORE | RUNNER_UP | NOTES |"
        genes = "marker_a;marker_b;marker_c;marker_d"
        assert header in output
        assert (
            f"| {DATA / 'exact.fa'} | pathotype | EHEC | {genes} | high"
            f" | 1.0000 | EPEC (0.0000) |  |"
        ) in output
        assert (
            f"| {DATA / 'exact.fa'} | toxin | - | {genes} | ambiguous | 1.0000 | - |"
            f" ambiguous: Toxin1 (1.0000), Toxin2 (1.0000) |"
        ) in output

    def test_v1_document_types_under_the_default_scheme_key(
        self, tmp_path_factory: pytest.TempPathFactory, tmp_path: Path
    ) -> None:
        """Given the same FASTA built with a typing/1 document, When typed,
        Then the single scheme evaluates identically and the phenotypes
        object keys on ``default``."""
        datadir = tmp_path_factory.mktemp("v1")
        result = runner.invoke(
            app,
            [
                "db",
                "build",
                DB,
                str(DATA / "markers.fa"),
                "--datadir",
                str(datadir),
                "--typing",
                str(DATA / "typing_markers_v1.json"),
            ],
        )
        assert result.exit_code == 0, result.stderr
        phenotypes = typing_json(datadir, "exact", tmp_path)["files"][0]["phenotypes"]
        assert list(phenotypes) == ["default"]
        assert phenotypes["default"]["phenotype"] == "EHEC"
        assert phenotypes["default"]["confidence"] == "high"


class TestBuildValidation:
    def test_unknown_gene_fails_the_build(self, tmp_path: Path) -> None:
        """Given a typing document referencing a gene the FASTA lacks, When
        built, Then exit 4 with TYPING_UNKNOWN_GENE and no database
        directory is left behind."""
        bad = tmp_path / "ghost.json"
        bad.write_text(
            json.dumps(
                {
                    "schema": "gapit.typing/2",
                    "schemes": [
                        {
                            "name": "pathotype",
                            "rules": [
                                {
                                    "model": "weighted_genes",
                                    "phenotype": "X",
                                    "weights": {"ghost": 1.0},
                                    "identity_floor": 90.0,
                                }
                            ],
                            "cutoff": 0.9,
                            "ambiguity_margin": 0.05,
                            "fallback": "unknown",
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )
        datadir = tmp_path / "datadir"
        datadir.mkdir()

        result = runner.invoke(
            app,
            [
                "db",
                "build",
                DB,
                str(DATA / "markers.fa"),
                "--datadir",
                str(datadir),
                "--typing",
                str(bad),
            ],
        )

        assert result.exit_code == 4
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "TYPING_UNKNOWN_GENE"
        assert envelope["context"]["gene"] == "ghost"
        assert not (datadir / DB).exists()

    def test_cluster_match_rule_is_rejected_on_a_gene_db(self, tmp_path: Path) -> None:
        """Given a typing document with a cluster_match rule (locus-level
        scoring), When built against a FASTA, Then exit 4 with
        TYPING_MALFORMED naming the rule model (DEC v1 is weighted_genes
        only)."""
        bad = tmp_path / "cm.json"
        bad.write_text(
            json.dumps(
                {
                    "schema": "gapit.typing/2",
                    "schemes": [
                        {
                            "name": "pathotype",
                            "rules": [
                                {
                                    "model": "cluster_match",
                                    "phenotype": "X",
                                    "coverage": {"weight": 0.5, "floor": 90.0},
                                    "identity": {"weight": 0.5, "floor": 90.0},
                                    "key_genes": {"weight": 0.0, "genes": ["marker_a"]},
                                }
                            ],
                            "cutoff": 0.9,
                            "ambiguity_margin": 0.05,
                            "fallback": "unknown",
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )
        datadir = tmp_path / "datadir"
        datadir.mkdir()

        result = runner.invoke(
            app,
            [
                "db",
                "build",
                DB,
                str(DATA / "markers.fa"),
                "--datadir",
                str(datadir),
                "--typing",
                str(bad),
            ],
        )

        assert result.exit_code == 4
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "TYPING_MALFORMED"
        assert envelope["context"]["model"] == "cluster_match"

    def test_manifest_records_the_typing_schema(self, typed_db: Path) -> None:
        """Given the typed build, When its manifest is read, Then
        typing_schema records the document's literal version."""
        manifest = json.loads((typed_db / DB / "gapit-manifest.json").read_text(encoding="utf-8"))
        assert manifest["typing_schema"] == "gapit.typing/2"
        assert manifest["kind"] == "gene"


class TestMcp:
    def test_mcp_screen_returns_the_pure_report(self, typed_db: Path) -> None:
        """Given the MCP screen tool and the typed gene db, When called,
        Then the gapit.report/1 document carries no phenotypes object —
        designation is not part of screening (a typing MCP tool is future
        work)."""
        raw = (
            json.dumps(
                {
                    "jsonrpc": "2.0",
                    "id": 1,
                    "method": "tools/call",
                    "params": {
                        "name": "screen",
                        "arguments": {
                            "files": [str(DATA / "exact.fa")],
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
        document = json.loads(response["result"]["content"][0]["text"])
        for entry in document["files"]:
            assert "phenotypes" not in entry


class TestGoldens:
    def test_typed_screen_goldens_at_pinned_now(
        self, typed_db: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Given the typed db and the exact sample, When the engine renders
        JSON and Markdown at the pinned now, Then each is byte-identical to
        the committed golden — a typed database screens the pure,
        phenotype-free report."""
        assert (typed_db / DB / "typing.json").is_file()
        ensure_blast()
        monkeypatch.chdir(DATA)
        report = screen_file(Path("exact.fa"), self._database(typed_db), PARAMS, dbtype="nucl")
        assert render_json([report], PARAMS, now=PINNED_NOW) == (
            GOLDEN / "gene_typing_exact.json"
        ).read_text(encoding="utf-8")
        assert render_markdown([report], PARAMS, now=PINNED_NOW) == (
            GOLDEN / "gene_typing_exact.md"
        ).read_text(encoding="utf-8")

    def test_typing_tsv_golden(
        self, typed_db: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Given the typed db and all four samples screened into one table
        (relative FILE keys via chdir), When typed as TSV, Then the output
        is byte-identical to the committed golden."""
        monkeypatch.chdir(DATA)
        table = tmp_path / "screen.tsv"
        result = runner.invoke(
            app,
            [
                "screen",
                *(f"{sample}.fa" for sample in SAMPLES),
                "--db",
                DB,
                "--datadir",
                str(typed_db),
                "-o",
                str(table),
                "--quiet",
            ],
        )
        assert result.exit_code == 0, result.stderr
        typed = runner.invoke(app, ["typing", str(table), "--datadir", str(typed_db), "--quiet"])
        assert typed.exit_code == 0, typed.stderr
        assert typed.stdout == (GOLDEN / "typing_markers.tsv").read_text(encoding="utf-8")

    @staticmethod
    def _database(datadir: Path):
        return next(db for db in discover_databases(datadir) if db.name == DB)
