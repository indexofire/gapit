"""Bundled lm_doumith: the Listeria Doumith-2004 serogrouping panel.

The sixth bundled database — five public-domain INSDC-extracted markers
(prs, lmo0737, lmo1118, ORF2819, ORF2110) behind the gapit.typing/2 scheme
``doumith_serogroup``: Doumith's multiplex-PCR table (Doumith et al. 2004,
J Clin Microbiol 42:3819) extended with the Huang 2011 4b variant (IVb-v,
declared first so the shared ORF2819+ORF2110 profile resolves to it). These
tests build the panel from the committed fixture into a tmp datadir (the
synthetic "isolates" are concatenations of marker subsets — zero network),
pin the designation matrix, the prs control-gene gate, the fallback paths,
and the two-stage golden, plus the bundled integration (auto-materialization
from the wheel copy, `db list` transition). Uppercase gene ids (ORF2819,
ORF2110) are part of the locked contract: build, screen, and typing carry
them verbatim.
"""

import json
from pathlib import Path

from typer.testing import CliRunner

from gapit.bundled import find_bundled
from gapit.cli import app
from gapit.fasta import iter_fasta
from gapit.records import read_manifest
from gapit.typing_models import read_typing_document
from gapit.typing_rules import ExactSetRule

SCHEMES = Path(__file__).parent / "data" / "typing" / "schemes"
GOLDEN = Path(__file__).parent / "golden"
PANEL = SCHEMES / "lm_doumith.fa"
SPEC = SCHEMES / "lm_doumith.json"
BUNDLE = Path(__file__).parents[1] / "src" / "gapit" / "data" / "dbs" / "lm_doumith"

runner = CliRunner()

FALLBACK = "untypeable (4a/4c, atypical profile, or non-Lm Listeria)"

# gene -> sequence from the committed panel (byte-equal to the wheel copy)
MARKERS = {record.id.split("~~~")[1]: record.sequence for record in iter_fasta(PANEL)}

# marker subset -> expected doumith_serogroup call (Doumith 2004 table 2;
# the all-five subset matches no published pattern and falls back)
CASES = [
    (["prs", "lmo0737"], "IIa"),
    (["prs", "lmo0737", "lmo1118"], "IIc"),
    (["prs", "ORF2819"], "IIb"),
    (["prs", "ORF2819", "ORF2110"], "IVb"),
    (["prs", "ORF2819", "ORF2110", "lmo0737"], "IVb-v"),
    (["prs"], FALLBACK),
    (["prs", "lmo0737", "lmo1118", "ORF2819", "ORF2110"], FALLBACK),
]


def build_db(datadir: Path) -> None:
    """`db build` the fixture panel + typing spec into a tmp datadir."""
    built = runner.invoke(
        app,
        [
            "db",
            "build",
            "lm_doumith",
            str(PANEL),
            "--typing",
            str(SPEC),
            "--datadir",
            str(datadir),
            "--quiet",
        ],
    )
    assert built.exit_code == 0, built.stderr


def write_sample(directory: Path, name: str, genes: list[str]) -> Path:
    """One synthetic isolate: the subset's markers concatenated into a contig."""
    path = directory / f"{name}.fa"
    path.write_text(f">{name}\n{''.join(MARKERS[gene] for gene in genes)}\n", encoding="utf-8")
    return path


def type_json(sample: Path, datadir: Path, tmp_path: Path) -> dict[str, object]:
    """The two-stage pipeline for one sample: screen -o table, then type it
    (gapit.typing_result/1); returns the scheme call (phenotype, notes, ...)."""
    table = tmp_path / f"{sample.stem}.tsv"
    screened = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            "lm_doumith",
            "--datadir",
            str(datadir),
            "-o",
            str(table),
            "--quiet",
        ],
    )
    assert screened.exit_code == 0, screened.stderr
    result = runner.invoke(
        app, ["typing", str(table), "--datadir", str(datadir), "--format", "json", "--quiet"]
    )
    assert result.exit_code == 0, result.stderr
    document = json.loads(result.stdout)
    assert document["db"] == "lm_doumith"
    return document["files"][0]["phenotypes"]["doumith_serogroup"]


def test_panel_and_spec_contract() -> None:
    """Given the committed panel and typing spec, When parsed, Then the five
    markers carry uppercase gene ids verbatim (ORF2819/ORF2110), the spec
    declares the one scheme in the pinned rule order with 95/95 floors and
    the prs control gate, and the wheel bundle is byte-equal to the fixture."""
    assert set(MARKERS) == {"prs", "lmo0737", "lmo1118", "ORF2819", "ORF2110"}
    scheme = read_typing_document(SPEC).schemes[0]
    assert scheme.name == "doumith_serogroup"
    assert [rule.phenotype for rule in scheme.rules] == ["IVb-v", "IVb", "IIb", "IIc", "IIa"]
    assert scheme.control_gene == "prs"
    assert scheme.cutoff == 1.0 and scheme.ambiguity_margin == 0.0
    assert scheme.fallback == FALLBACK
    for rule in scheme.rules:
        assert isinstance(rule, ExactSetRule)
        assert rule.identity_floor == 95.0 and rule.coverage_floor == 95.0
    assert (BUNDLE / "sequences").read_bytes() == PANEL.read_bytes()
    assert (BUNDLE / "typing.json").read_bytes() == SPEC.read_bytes()


def test_doumith_serogroup_matrix(tmp_path: Path) -> None:
    """Given the built panel database and synthetic marker-subset isolates,
    When each runs the two-stage pipeline, Then every Doumith call holds —
    the four serogroup patterns, the IVb-v variant, and the two fallbacks
    (prs only; all five markers = no published pattern)."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    build_db(datadir)
    for genes, expected in CASES:
        sample = write_sample(tmp_path, "lm_" + "_".join(genes), genes)
        call = type_json(sample, datadir, tmp_path)
        assert call["phenotype"] == expected, genes
        assert call["confidence"] == ("high" if expected != FALLBACK else "low")


def test_prs_control_gate_and_unrelated_sequence(tmp_path: Path) -> None:
    """Given markers without prs and an unrelated sequence, When typed, Then
    the prs-less profile falls back with the ``control gene absent`` note
    (prs is genus-level, not Lm-specific), and the unrelated sequence yields
    zero screen rows, so `typing` refuses with TYPING_NO_DATA (exit 5)."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    build_db(datadir)

    no_prs = write_sample(tmp_path, "lm_no_prs", ["lmo0737", "ORF2819"])
    call = type_json(no_prs, datadir, tmp_path)
    assert call["phenotype"] == FALLBACK
    assert call["notes"] == ["control gene absent"]

    unrelated = tmp_path / "lm_unrelated.fa"
    unrelated.write_text(
        ">lm_unrelated\n" + "ACGTTGCAAGGCTTACGGATCCTTAGGCATCGGAATTCGGCTAACGGGATCC" * 20 + "\n",
        encoding="utf-8",
    )
    table = tmp_path / "unrelated.tsv"
    screened = runner.invoke(
        app,
        [
            "screen",
            str(unrelated),
            "--db",
            "lm_doumith",
            "--datadir",
            str(datadir),
            "-o",
            str(table),
            "--quiet",
        ],
    )
    assert screened.exit_code == 0, screened.stderr
    result = runner.invoke(
        app, ["typing", str(table), "--datadir", str(datadir), "--format", "json", "--quiet"]
    )
    assert result.exit_code == 5
    assert json.loads(result.stderr)["code"] == "TYPING_NO_DATA"


def test_typing_command_end_to_end_on_one_representative(tmp_path: Path) -> None:
    """Given one representative IIa isolate, When `screen -o table` writes
    the gene table and `gapit typing <table>` designates from that file arg,
    Then the default TSV renders the eight columns with the IIa call and the
    file's sorted uppercase-carrying gene list."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    build_db(datadir)
    sample = write_sample(tmp_path, "lm_doumith_iia", ["prs", "lmo0737"])
    table = tmp_path / "iia.tsv"
    screened = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            "lm_doumith",
            "--datadir",
            str(datadir),
            "-o",
            str(table),
            "--nopath",
            "--quiet",
        ],
    )
    assert screened.exit_code == 0, screened.stderr
    result = runner.invoke(app, ["typing", str(table), "--datadir", str(datadir), "--quiet"])
    assert result.exit_code == 0, result.stderr
    lines = result.stdout.splitlines()
    assert lines[0].split("\t") == [
        "FILE",
        "SCHEME",
        "PHENOTYPE",
        "GENES",
        "CONFIDENCE",
        "SCORE",
        "RUNNER_UP",
        "NOTES",
    ]
    row = lines[1].split("\t")
    assert row[:4] == ["lm_doumith_iia.fa", "doumith_serogroup", "IIa", "lmo0737;prs"]


def test_typing_golden_iia_and_ivb(tmp_path: Path) -> None:
    """Given the IIa and IVb isolates screened into one table, When typed as
    the DEFAULT TSV, Then the designation output is byte-identical to the
    committed golden — the two-stage lm_doumith surface, frozen."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    build_db(datadir)
    iia = write_sample(tmp_path, "lm_doumith_iia", ["prs", "lmo0737"])
    ivb = write_sample(tmp_path, "lm_doumith_ivb", ["prs", "ORF2819", "ORF2110"])
    table = tmp_path / "golden.tsv"
    screened = runner.invoke(
        app,
        [
            "screen",
            str(iia),
            str(ivb),
            "--db",
            "lm_doumith",
            "--datadir",
            str(datadir),
            "--nopath",
            "-o",
            str(table),
            "--quiet",
        ],
    )
    assert screened.exit_code == 0, screened.stderr
    result = runner.invoke(app, ["typing", str(table), "--datadir", str(datadir), "--quiet"])
    assert result.exit_code == 0, result.stderr
    assert result.stdout == (GOLDEN / "lm_doumith_typing.tsv").read_text(encoding="utf-8")


def test_bundled_lm_doumith_auto_materializes_and_lists(tmp_path: Path) -> None:
    """Given a fresh EMPTY datadir and no prior fetch, When screening against
    the bundled lm_doumith, Then the wheel copy materializes (5 records, one
    stderr note), the manifest certifies the bundled-source typed database,
    every marker is called at 100%, and `db list` flips the row from bundled
    to installed (5)."""
    assert find_bundled("lm_doumith") is not None
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    sample = write_sample(tmp_path, "lm_doumith_ivbv", ["prs", "ORF2819", "ORF2110", "lmo0737"])

    result = runner.invoke(
        app, ["screen", str(sample), "--db", "lm_doumith", "--datadir", str(datadir), "--quiet"]
    )
    assert result.exit_code == 0, result.stderr
    manifest = read_manifest(datadir / "lm_doumith" / "gapit-manifest.json")
    assert manifest.n_records == 5
    assert manifest.source == "bundled"
    assert manifest.typing_schema == "gapit.typing/2"
    assert (datadir / "lm_doumith" / "typing.json").is_file()
    assert (datadir / "lm_doumith" / "sequences.nin").is_file()

    listing = runner.invoke(app, ["db", "list", "--datadir", str(datadir)])
    assert listing.exit_code == 0
    assert (
        "lm_doumith\tgapit-curated (public-domain INSDC sources)\tinstalled (5)\tnucl\t"
        "Listeria monocytogenes serogrouping (Doumith 2004)" in listing.stdout
    )

    table = tmp_path / "bundled.tsv"
    assert (
        runner.invoke(
            app,
            [
                "screen",
                str(sample),
                "--db",
                "lm_doumith",
                "--datadir",
                str(datadir),
                "-o",
                str(table),
                "--nopath",
                "--quiet",
            ],
        ).exit_code
        == 0
    )
    typed = runner.invoke(
        app, ["typing", str(table), "--datadir", str(datadir), "--format", "json", "--quiet"]
    )
    assert typed.exit_code == 0, typed.stderr
    call = json.loads(typed.stdout)["files"][0]["phenotypes"]["doumith_serogroup"]
    assert call["phenotype"] == "IVb-v"
