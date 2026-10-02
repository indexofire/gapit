"""Gene-path typing stage-2 integration: the rule primitives end to end.

A typed marker fixture database (``serotyping``: Doumith-style Listeria
exact_set table gated by the prs control gene, a two-wzx o_group scheme
with the unique-group mixed override, a weighted k_group scheme, and the
``{o_group}:{k_group}`` compose scheme) built with ``db build --typing``
and run through the two-stage pipeline: ``screen -o result.tsv`` (real
blastn engine), then ``gapit typing result.tsv``.

Covers: the exact_set Doumith calls with surfaced notes, control-gate
zeroing, mixed-infection output with the pair in ambiguous, compose with
both ingredients called versus one ambiguous, the typed/untyped screen
byte-identity (the purification pin), build-time scheme-gene validation,
and goldens at a pinned timestamp. Unit math lives in
test_typing_stage2.py; schema validation in test_typing_models.py; the
command's input parsing in test_typing_command.py.
"""

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from gapit.blast import ensure_blast, screen_file
from gapit.cli import app
from gapit.db import discover_databases
from gapit.formats.json import render_json
from gapit.formats.md import render_markdown
from gapit.report import ScreeningParams

DATA = Path(__file__).parent / "data" / "typing"
GOLDEN = Path(__file__).parent / "golden"
DB = "serotyping"
PINNED_NOW = datetime(2026, 10, 1, 12, 0, 0, tzinfo=UTC)
PARAMS = ScreeningParams(db=DB)

runner = CliRunner()

SAMPLES = ("listeria_4b", "listeria_1c", "listeria_nogate", "vp_typed", "vp_mixed")


@pytest.fixture(scope="module")
def typed_db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The serotyping fixture db with the four-scheme typing/2 document,
    plus an untyped rebuild of the same FASTA in a sibling datadir."""
    datadir = tmp_path_factory.mktemp("serotyped")
    result = runner.invoke(
        app,
        [
            "db",
            "build",
            DB,
            str(DATA / "serotyping.fa"),
            "--datadir",
            str(datadir),
            "--typing",
            str(DATA / "typing_serotyping_v2.json"),
        ],
    )
    assert result.exit_code == 0, result.stderr
    return datadir


@pytest.fixture(scope="module")
def untyped_db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    datadir = tmp_path_factory.mktemp("serotyped_untyped")
    result = runner.invoke(
        app, ["db", "build", DB, str(DATA / "serotyping.fa"), "--datadir", str(datadir)]
    )
    assert result.exit_code == 0, result.stderr
    return datadir


def typing_json(datadir: Path, sample: str, tmp_path: Path) -> dict[str, Any]:
    """The typing_result/1 document for one sample: screen it to a table
    (stage 1), then type the table (stage 2)."""
    table = tmp_path / f"{sample}.tsv"
    screened = runner.invoke(
        app,
        [
            "screen",
            str(DATA / f"{sample}.fa"),
            "--db",
            DB,
            "--datadir",
            str(datadir),
            "-o",
            str(table),
            "--quiet",
        ],
    )
    assert screened.exit_code == 0, screened.stderr
    result = runner.invoke(
        app,
        [
            "typing",
            str(table),
            "--datadir",
            str(datadir),
            "--format",
            "json",
            "--quiet",
        ],
    )
    assert result.exit_code == 0, result.stderr
    return json.loads(result.stdout)


class TestDoumithTable:
    def test_4b_sample_calls_the_exact_set_rule_with_notes(
        self, typed_db: Path, tmp_path: Path
    ) -> None:
        """Given the 4b sample (prs + orf2110 + orf2870 at 100%), When
        typed, Then doumith calls 4b-4d-4e at 1.0 with the rule's notes
        verbatim and one 1.0 component per required gene."""
        phenotypes = typing_json(typed_db, "listeria_4b", tmp_path)["files"][0]["phenotypes"]
        doumith = phenotypes["doumith"]
        assert doumith["phenotype"] == "4b-4d-4e"
        assert doumith["confidence"] == "high"
        assert doumith["score"] == 1.0
        assert doumith["components"] == [
            {"name": "orf2110", "score": 1.0},
            {"name": "orf2870", "score": 1.0},
        ]
        assert doumith["notes"] == ["4b/4d/4e share one PCR pattern (Doumith et al. 2004)"]

    def test_1c_sample_calls_the_excluded_marker_rule(self, typed_db: Path, tmp_path: Path) -> None:
        """Given the 1c sample (prs + lmo1118 + lmo0733), When typed, Then
        doumith calls 1/2c-3c — the excluded lmo0737 is absent and both
        required markers are present."""
        phenotypes = typing_json(typed_db, "listeria_1c", tmp_path)["files"][0]["phenotypes"]
        doumith = phenotypes["doumith"]
        assert doumith["phenotype"] == "1/2c-3c"
        assert doumith["confidence"] == "high"
        assert doumith["components"] == [
            {"name": "lmo1118", "score": 1.0},
            {"name": "lmo0733", "score": 1.0},
            {"name": "lmo0737", "score": 1.0},
        ]

    def test_fallback_when_no_rule_matches(self, typed_db: Path, tmp_path: Path) -> None:
        """Given the typed vp sample (no Listeria markers except the panel
        genes it carries), When typed, Then doumith falls back to NT — via
        the control gate (prs absent) with the control note."""
        phenotypes = typing_json(typed_db, "vp_typed", tmp_path)["files"][0]["phenotypes"]
        doumith = phenotypes["doumith"]
        assert doumith["phenotype"] == "NT"
        assert doumith["confidence"] == "low"
        assert doumith["notes"] == ["control gene absent"]


class TestControlGate:
    def test_missing_control_gene_zeroes_the_scheme(self, typed_db: Path, tmp_path: Path) -> None:
        """Given the nogate sample (lmo0737 present but prs absent), When
        typed, Then doumith outputs NT with low confidence, score 0.0, and
        the 'control gene absent' note — the 1/2a-3a markers alone cannot
        call."""
        phenotypes = typing_json(typed_db, "listeria_nogate", tmp_path)["files"][0]["phenotypes"]
        doumith = phenotypes["doumith"]
        assert doumith["phenotype"] == "NT"
        assert doumith["confidence"] == "low"
        assert doumith["score"] == 0.0
        assert doumith["components"] == []
        assert doumith["notes"] == ["control gene absent"]


class TestMixedInfection:
    def test_two_wzx_genes_call_mixed(self, typed_db: Path, tmp_path: Path) -> None:
        """Given the mixed sample (wzx_o1 AND wzx_o2 present), When typed,
        Then o_group outputs the mixed phenotype with the pair listed in
        ambiguous and the mixed-infection note."""
        phenotypes = typing_json(typed_db, "vp_mixed", tmp_path)["files"][0]["phenotypes"]
        o_group = phenotypes["o_group"]
        assert o_group["phenotype"] == "mixed"
        assert o_group["confidence"] == "low"
        assert o_group["ambiguous"] == [
            {"phenotype": "wzx_o1", "score": 1.0},
            {"phenotype": "wzx_o2", "score": 1.0},
        ]
        assert o_group["notes"] == ["mixed infection in unique group 'wzx'"]

    def test_single_wzx_gene_calls_the_o_group(self, typed_db: Path, tmp_path: Path) -> None:
        """Given the typed sample (wzx_o1 only), When typed, Then o_group
        calls O1 with high confidence and no mixed override."""
        phenotypes = typing_json(typed_db, "vp_typed", tmp_path)["files"][0]["phenotypes"]
        o_group = phenotypes["o_group"]
        assert o_group["phenotype"] == "O1"
        assert o_group["confidence"] == "high"
        assert "ambiguous" not in o_group


class TestCompose:
    def test_both_ingredients_called_compose_high(self, typed_db: Path, tmp_path: Path) -> None:
        """Given O1 and K13 both called with high confidence, When typed,
        Then serotype renders O1:K13 at the minimum ingredient score with
        one component per ingredient."""
        phenotypes = typing_json(typed_db, "vp_typed", tmp_path)["files"][0]["phenotypes"]
        serotype = phenotypes["serotype"]
        assert serotype["phenotype"] == "O1:K13"
        assert serotype["confidence"] == "high"
        assert serotype["score"] == 1.0
        assert serotype["components"] == [
            {"name": "o_group", "score": 1.0},
            {"name": "k_group", "score": 1.0},
        ]

    def test_fallback_ingredients_flow_through(self, typed_db: Path, tmp_path: Path) -> None:
        """Given a Listeria sample where o_group/k_group fall back to their
        strings, When typed, Then serotype renders O?:K? with low confidence
        (the fallback strings flow through the template)."""
        phenotypes = typing_json(typed_db, "listeria_4b", tmp_path)["files"][0]["phenotypes"]
        serotype = phenotypes["serotype"]
        assert serotype["phenotype"] == "O?:K?"
        assert serotype["confidence"] == "low"
        assert serotype["score"] == 0.0

    def test_ambiguous_ingredient_makes_the_composition_null(
        self, typed_db: Path, tmp_path: Path
    ) -> None:
        """Given the mixed sample (k_group ambiguous on the K13/K64 tie),
        When typed, Then serotype is null with k_group's variants in
        ambiguous and a note naming the ingredient."""
        phenotypes = typing_json(typed_db, "vp_mixed", tmp_path)["files"][0]["phenotypes"]
        assert phenotypes["k_group"]["phenotype"] is None
        serotype = phenotypes["serotype"]
        assert serotype["phenotype"] is None
        assert serotype["confidence"] == "ambiguous"
        assert serotype["ambiguous"] == [
            {"phenotype": "K13", "score": 1.0},
            {"phenotype": "K64", "score": 1.0},
        ]
        assert serotype["notes"] == ["ingredient scheme 'k_group' is ambiguous"]


class TestPurifiedScreen:
    def test_typed_and_untyped_screens_are_byte_identical(
        self, typed_db: Path, untyped_db: Path
    ) -> None:
        """Given the same inputs against the typed and the untyped build,
        When screened in the default format, Then the outputs are
        byte-identical 15-column tables — designation lives in `gapit
        typing` (the inverted former phenotype-column pin)."""
        inputs = [str(DATA / f"{sample}.fa") for sample in SAMPLES]
        typed = runner.invoke(
            app, ["screen", *inputs, "--db", DB, "--datadir", str(typed_db), "--quiet"]
        )
        untyped = runner.invoke(
            app, ["screen", *inputs, "--db", DB, "--datadir", str(untyped_db), "--quiet"]
        )
        assert typed.exit_code == 0, typed.stderr
        assert untyped.exit_code == 0, untyped.stderr
        assert typed.stdout == untyped.stdout
        for line in untyped.stdout.splitlines():
            assert len(line.split("\t")) == 15

    def test_untyped_json_has_no_phenotypes(self, untyped_db: Path) -> None:
        """Given the untyped rebuild, When screened as JSON, Then no file
        block carries a phenotypes object."""
        for sample in ("listeria_4b", "vp_mixed"):
            result = runner.invoke(
                app,
                [
                    "screen",
                    str(DATA / f"{sample}.fa"),
                    "--db",
                    DB,
                    "--datadir",
                    str(untyped_db),
                    "--format",
                    "json",
                    "--quiet",
                ],
            )
            assert result.exit_code == 0, result.stderr
            for entry in json.loads(result.stdout)["files"]:
                assert "phenotypes" not in entry


class TestBuildValidation:
    def test_unknown_control_gene_fails_the_build(self, tmp_path: Path) -> None:
        """Given a typing document whose scheme control gene the FASTA
        lacks, When built, Then exit 4 with TYPING_UNKNOWN_GENE naming the
        scheme and no database directory is left behind."""
        bad = tmp_path / "bad.json"
        bad.write_text(
            json.dumps(
                {
                    "schema": "gapit.typing/2",
                    "schemes": [
                        {
                            "name": "doumith",
                            "rules": [
                                {
                                    "model": "exact_set",
                                    "phenotype": "X",
                                    "requires": ["prs"],
                                }
                            ],
                            "cutoff": 0.9,
                            "ambiguity_margin": 0.05,
                            "fallback": "NT",
                            "control_gene": "ipaH",
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
                str(DATA / "serotyping.fa"),
                "--datadir",
                str(datadir),
                "--typing",
                str(bad),
            ],
        )

        assert result.exit_code == 4
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "TYPING_UNKNOWN_GENE"
        assert envelope["context"]["scheme"] == "doumith"
        assert not (datadir / DB).exists()


class TestGoldens:
    def test_goldens_at_pinned_now(self, typed_db: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Given the typed db and the 4b/mixed samples, When the engine
        renders the screen report as JSON and Markdown at the pinned now,
        Then each output is byte-identical to its committed golden — the
        typed database screens the pure, phenotype-free report."""
        assert (typed_db / DB / "typing.json").is_file()
        ensure_blast()
        monkeypatch.chdir(DATA)
        for sample, stem in (
            ("listeria_4b", "gene_serotyping_l4b"),
            ("vp_mixed", "gene_serotyping_mixed"),
        ):
            report = screen_file(
                Path(f"{sample}.fa"), self._database(typed_db), PARAMS, dbtype="nucl"
            )
            assert render_json([report], PARAMS, now=PINNED_NOW) == (
                GOLDEN / f"{stem}.json"
            ).read_text(encoding="utf-8")
            assert render_markdown([report], PARAMS, now=PINNED_NOW) == (
                GOLDEN / f"{stem}.md"
            ).read_text(encoding="utf-8")

    def test_typing_tsv_golden(
        self, typed_db: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Given the typed db and all five samples screened into one table
        (relative FILE keys via chdir), When typed as TSV, Then the output
        is byte-identical to the committed golden — compose, mixed, and the
        control-gate fallback in one seven-column table."""
        monkeypatch.chdir(DATA)
        table = tmp_path / "screen.tsv"
        screened = runner.invoke(
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
        assert screened.exit_code == 0, screened.stderr
        result = runner.invoke(app, ["typing", str(table), "--datadir", str(typed_db), "--quiet"])
        assert result.exit_code == 0, result.stderr
        assert result.stdout == (GOLDEN / "typing_serotyping.tsv").read_text(encoding="utf-8")

    @staticmethod
    def _database(datadir: Path):
        return next(db for db in discover_databases(datadir) if db.name == DB)
