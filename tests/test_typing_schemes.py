"""typing/2 stage 3 — the gene-path scheme cookbook, self-checked.

The researched designation schemes as validated example documents
(``tests/data/typing/schemes/``): Doumith's complete Listeria table, the
ShigaTyper-semantics Shigella/EIEC skeleton, the meningotype serogroup
panel with the EX7E allele-probe trick, and the dual-scheme DEC
designation (GB 4789.6-2016 vs the risk-monitoring variant). Each fixture
db is built with ``db build --typing`` and run through the two-stage
pipeline — ``screen -o table.tsv`` (real blastn engine), then ``gapit
typing table.tsv`` — with samples deterministic concatenations/mutations
of the fixture markers (synthetic content only — curation is future
work). These are the equivalence locks: the typing command's calls equal
the calls the inline engine produced before designation moved to the
command. Document round-trip validation covers all six cookbook files
(the two cluster ones run in test_typing_schemes_cluster.py).
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
from gapit.fasta import iter_fasta
from gapit.formats.json import render_json
from gapit.formats.md import render_markdown
from gapit.report import ScreeningParams
from gapit.typing_models import read_typing_document, template_placeholders

SCHEMES = Path(__file__).parent / "data" / "typing" / "schemes"
TYPING = Path(__file__).parent / "data" / "typing"
GOLDEN = Path(__file__).parent / "golden"
PINNED_NOW = datetime(2026, 10, 1, 12, 0, 0, tzinfo=UTC)
FLIP = {"A": "C", "C": "G", "G": "T", "T": "A"}

runner = CliRunner()

GENE_FIXTURES = ("doumith", "shigella", "meningotype", "dec")
COOKBOOK_DOCS = ("doumith", "shigella", "meningotype", "dec", "vp_ok", "cholerae")


def mutate(parent: str, *sites: int) -> str:
    bases = list(parent)
    for site in sites:
        bases[site] = FLIP[bases[site]]
    return "".join(bases)


def write_sample(directory: Path, name: str, contigs: list[tuple[str, str]]) -> Path:
    path = directory / f"{name}.fa"
    path.write_text(
        "".join(f">{contig}\n{sequence}\n" for contig, sequence in contigs), encoding="utf-8"
    )
    return path


def markers(fasta: Path) -> dict[str, str]:
    return {record.id: record.sequence for record in iter_fasta(fasta)}


def screen_json(datadir: Path, db: str, sample: Path) -> dict[str, Any]:
    result = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            db,
            "--datadir",
            str(datadir),
            "--format",
            "json",
            "--quiet",
        ],
    )
    assert result.exit_code == 0, result.stderr
    return json.loads(result.stdout)


def typing_document(datadir: Path, db: str, sample: Path, tmp_path: Path) -> dict[str, Any]:
    """The typing_result/1 document for one sample — the two-stage pipeline
    (screen to a table, then type it) end to end through the CLI."""
    table = tmp_path / f"{sample.stem}.tsv"
    screened = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            db,
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
        ["typing", str(table), "--datadir", str(datadir), "--format", "json", "--quiet"],
    )
    assert result.exit_code == 0, result.stderr
    return json.loads(result.stdout)


def call(document: dict[str, Any], scheme: str) -> dict[str, Any]:
    return document["files"][0]["phenotypes"][scheme]


@pytest.fixture(scope="module")
def gene_dbs(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The four gene-kind cookbook databases built into one fresh datadir
    (each with its full typing/2 document)."""
    datadir = tmp_path_factory.mktemp("schemes_gene")
    for name in GENE_FIXTURES:
        result = runner.invoke(
            app,
            [
                "db",
                "build",
                name,
                str(SCHEMES / f"{name}.fa"),
                "--datadir",
                str(datadir),
                "--typing",
                str(SCHEMES / f"{name}.json"),
            ],
        )
        assert result.exit_code == 0, result.stderr
    return datadir


class TestDocumentsValidate:
    @pytest.mark.parametrize("name", COOKBOOK_DOCS)
    def test_round_trip_load(self, name: str) -> None:
        """Given any cookbook typing document, When read through the typing
        loader, Then it validates against the current gapit.typing/2 schema
        (the committed fixtures can never drift from the model)."""
        path = SCHEMES / f"{name}.json"
        document = read_typing_document(path)
        assert document.schema_name == "gapit.typing/2"
        assert document.schemes

    def test_vp_compose_placeholders_name_the_sibling_schemes(self) -> None:
        """Given the vp_ok document, When its compose template is parsed,
        Then the placeholders are exactly the o_group/k_group siblings in
        order."""
        document = read_typing_document(SCHEMES / "vp_ok.json")
        (compose,) = [scheme.compose for scheme in document.schemes if scheme.compose]
        assert template_placeholders(compose) == ("o_group", "k_group")


class TestDoumith:
    @staticmethod
    def samples(tmp_path: Path) -> dict[str, Path]:
        g = markers(SCHEMES / "doumith.fa")
        return {
            name: write_sample(tmp_path, name, [(name, "".join(g[marker] for marker in panel))])
            for name, panel in {
                "lm_1_2a": ["prs", "lmo0737"],
                "lm_1_2c": ["prs", "lmo0737", "lmo1118"],
                "lm_1_2b": ["prs", "orf2819"],
                "lm_4b": ["prs", "orf2819", "orf2110"],
                "lm_4b_star": ["prs", "lmo0737", "orf2819", "orf2110"],
            }.items()
        }

    @pytest.mark.parametrize(
        ("sample", "expected"),
        [
            ("lm_1_2a", "1/2a-3a"),
            ("lm_1_2c", "1/2c-3c"),
            ("lm_1_2b", "1/2b-3b-7"),
            ("lm_4b", "4b-4d-4e"),
        ],
    )
    def test_serogroup_profiles_call_their_rule(
        self, gene_dbs: Path, tmp_path: Path, sample: str, expected: str
    ) -> None:
        """Given a synthetic genome carrying exactly one Doumith marker
        profile, When screened, Then the doumith scheme calls that
        serogroup with high confidence at 1.0."""
        document = typing_document(gene_dbs, "doumith", self.samples(tmp_path)[sample], tmp_path)
        assert call(document, "doumith")["phenotype"] == expected

    def test_4b_call_carries_the_hgt_caveat_notes(self, gene_dbs: Path, tmp_path: Path) -> None:
        """Given the 4b sample, When screened, Then the winning rule's notes
        surface verbatim — the shared-pattern note plus the 4b* HGT caveat."""
        document = typing_document(gene_dbs, "doumith", self.samples(tmp_path)["lm_4b"], tmp_path)
        detail = call(document, "doumith")
        assert detail["score"] == 1.0
        assert detail["components"] == [
            {"name": "orf2819", "score": 1.0},
            {"name": "orf2110", "score": 1.0},
            {"name": "lmo0737", "score": 1.0},
            {"name": "lmo1118", "score": 1.0},
        ]
        assert "4b* (IVb-v1)" in detail["notes"][1]

    def test_4b_star_profile_falls_back_to_nt(self, gene_dbs: Path, tmp_path: Path) -> None:
        """Given the 4b* profile (every marker, the HGT-acquired lmo0737
        cassette), When screened, Then no plain pattern matches and the
        scheme falls back to NT with low confidence (the caveat working)."""
        document = typing_document(
            gene_dbs, "doumith", self.samples(tmp_path)["lm_4b_star"], tmp_path
        )
        detail = call(document, "doumith")
        assert detail["phenotype"] == "NT"
        assert detail["confidence"] == "low"

    def test_missing_prs_control_gates_the_scheme(self, gene_dbs: Path, tmp_path: Path) -> None:
        """Given a genome with lmo0737 but no prs, When screened, Then the
        control gate zeroes the scheme to NT with the control note."""
        g = markers(SCHEMES / "doumith.fa")
        sample = write_sample(tmp_path, "lm_nogate", [("lm_nogate", g["lmo0737"])])
        document = typing_document(gene_dbs, "doumith", sample, tmp_path)
        detail = call(document, "doumith")
        assert detail["phenotype"] == "NT"
        assert detail["notes"] == ["control gene absent"]


class TestShigella:
    @staticmethod
    def samples(tmp_path: Path) -> dict[str, Path]:
        g = markers(SCHEMES / "shigella.fa")
        return {
            name: write_sample(tmp_path, name, [(name, "".join(g[marker] for marker in panel))])
            for name, panel in {
                "sh_form_i": ["ipaH_c", "ipaB", "Ss_wzx", "Ss_wzy", "Ss_methylase"],
                "sh_form_ii": ["ipaH_c", "Ss_methylase"],
                "sh_dys1": ["ipaH_c", "Sd1_wzx", "Sd1_rfp"],
                "sh_flex_2a": ["ipaH_c", "Sf_wzx", "Sf_wzy", "gtrII"],
                "sh_flex_y": ["ipaH_c", "Sf_wzx", "Sf_wzy"],
                "sh_eiec": ["ipaH_c", "EclacY"],
                "sh_eiec_exempt": ["ipaH_c", "EclacY", "Sb9_wzx"],
                "sh_mixed": ["ipaH_c", "Ss_wzx", "Ss_wzy", "Sf_wzx", "Sf_wzy"],
            }.items()
        }

    def test_form_i_with_methylase_wins_by_declaration_order(
        self, gene_dbs: Path, tmp_path: Path
    ) -> None:
        """Given a form I sample that also carries Ss_methylase (the
        phase-variable antigen OFF), When screened, Then both the form I
        exact_set and the form II weighted rule score 1.0 and the declared
        form I rule wins the tie (the exact_set exception)."""
        document = typing_document(
            gene_dbs, "shigella", self.samples(tmp_path)["sh_form_i"], tmp_path
        )
        detail = call(document, "shigella")
        assert detail["phenotype"] == "S. sonnei form I"
        assert detail["confidence"] == "high"

    def test_form_ii_calls_through_the_methylase_alone(
        self, gene_dbs: Path, tmp_path: Path
    ) -> None:
        """Given ipaH_c + Ss_methylase and no wzx at all, When screened,
        Then the scheme calls form II."""
        document = typing_document(
            gene_dbs, "shigella", self.samples(tmp_path)["sh_form_ii"], tmp_path
        )
        assert call(document, "shigella")["phenotype"] == "S. sonnei form II"

    @pytest.mark.parametrize(
        ("sample", "expected"),
        [
            ("sh_dys1", "S. dysenteriae 1"),
            ("sh_flex_2a", "S. flexneri 2a"),
            ("sh_flex_y", "S. flexneri Y/novel"),
        ],
    )
    def test_serotype_rules_call(
        self, gene_dbs: Path, tmp_path: Path, sample: str, expected: str
    ) -> None:
        """Given the dysenteriae-1 and flexneri profiles, When screened,
        Then their exact_set rules call (Y/novel = the Sf base with no
        listed conversion gene)."""
        document = typing_document(gene_dbs, "shigella", self.samples(tmp_path)[sample], tmp_path)
        assert call(document, "shigella")["phenotype"] == expected

    def test_eiec_lacy_call_with_the_approximation_note(
        self, gene_dbs: Path, tmp_path: Path
    ) -> None:
        """Given ipaH_c + EclacY, When screened, Then the scheme calls EIEC
        and surfaces the approximation note (lacY+ minus the boydii 9/15
        exemptions)."""
        document = typing_document(
            gene_dbs, "shigella", self.samples(tmp_path)["sh_eiec"], tmp_path
        )
        detail = call(document, "shigella")
        assert detail["phenotype"] == "EIEC"
        assert "Approximation of ShigaTyper checkpoint 3" in detail["notes"][0]

    def test_exempt_boydii_9_blocks_the_eiec_call(self, gene_dbs: Path, tmp_path: Path) -> None:
        """Given lacY+ with Sb9_wzx present (the exemption), When screened,
        Then the EIEC exact_set fails on its exclude and the scheme falls
        back rather than calling EIEC."""
        document = typing_document(
            gene_dbs, "shigella", self.samples(tmp_path)["sh_eiec_exempt"], tmp_path
        )
        detail = call(document, "shigella")
        assert detail["phenotype"] != "EIEC"

    def test_two_wzx_call_mixed(self, gene_dbs: Path, tmp_path: Path) -> None:
        """Given a sample with both Ss_wzx and Sf_wzx, When screened, Then
        the wzx unique group fires the mixed phenotype with the pair in
        ambiguous and the mixed-infection note."""
        document = typing_document(
            gene_dbs, "shigella", self.samples(tmp_path)["sh_mixed"], tmp_path
        )
        detail = call(document, "shigella")
        assert detail["phenotype"] == "mixed Shigella serotypes"
        assert detail["confidence"] == "low"
        assert detail["ambiguous"] == [
            {"phenotype": "Ss_wzx", "score": 1.0},
            {"phenotype": "Sf_wzx", "score": 1.0},
        ]

    def test_absent_ipah_control_gates_the_scheme(self, gene_dbs: Path, tmp_path: Path) -> None:
        """Given EclacY without the ipaH_c control, When screened, Then the
        scheme outputs the not-Shigella-or-EIEC fallback with the control
        note (an E. coli that is not EIEC)."""
        g = markers(SCHEMES / "shigella.fa")
        sample = write_sample(tmp_path, "sh_negative", [("sh_negative", g["EclacY"])])
        document = typing_document(gene_dbs, "shigella", sample, tmp_path)
        detail = call(document, "shigella")
        assert detail["phenotype"] == "not Shigella or EIEC"
        assert detail["notes"] == ["control gene absent"]


class TestMeningotype:
    @staticmethod
    def samples(tmp_path: Path) -> dict[str, Path]:
        g = markers(SCHEMES / "meningotype.fa")
        # dual W/Y (EX7E = S): one SNP shared with each allele probe
        s_allele = mutate(g["synG_EX7E_P"], 47, 93)
        return {
            "nm_b": write_sample(tmp_path, "nm_b", [("nm_b", g["ctrA"] + g["synD"])]),
            "nm_w": write_sample(tmp_path, "nm_w", [("nm_w", g["ctrA"] + g["synG_EX7E_P"])]),
            "nm_y": write_sample(
                tmp_path,
                "nm_y",
                [("nm_y", g["ctrA"] + g["synF"] + g["synG_EX7E_G"])],
            ),
            "nm_dual": write_sample(tmp_path, "nm_dual", [("nm_dual", g["ctrA"] + s_allele)]),
            "nm_ng": write_sample(tmp_path, "nm_ng", [("nm_ng", g["ctrA"])]),
        }

    @pytest.mark.parametrize(("sample", "expected"), [("nm_b", "B"), ("nm_w", "W"), ("nm_y", "Y")])
    def test_serogroup_panel_calls(
        self, gene_dbs: Path, tmp_path: Path, sample: str, expected: str
    ) -> None:
        """Given synthetic serogroup profiles, When screened, Then the
        panel calls the serogroup — W and Y through their EX7E allele
        probes at the 99.5 identity floor."""
        document = typing_document(
            gene_dbs, "meningotype", self.samples(tmp_path)[sample], tmp_path
        )
        assert call(document, "serogroup")["phenotype"] == expected

    def test_dual_wy_strain_falls_back_to_ng(self, gene_dbs: Path, tmp_path: Path) -> None:
        """Given the dual W/Y S allele (one SNP from each probe), When
        screened, Then neither probe clears the 99.5 floor and the scheme
        falls back to NG — the allele-level ambiguity an allele_match
        primitive would resolve."""
        document = typing_document(
            gene_dbs, "meningotype", self.samples(tmp_path)["nm_dual"], tmp_path
        )
        detail = call(document, "serogroup")
        assert detail["phenotype"] == "NG"
        assert detail["confidence"] == "low"

    def test_ctrA_only_is_ng_and_no_ctrA_gates(self, gene_dbs: Path, tmp_path: Path) -> None:
        """Given ctrA alone, When screened, Then the scheme outputs NG below
        cutoff; given a non-meningococcal genome (no ctrA), the control gate
        adds the control note."""
        ng = typing_document(gene_dbs, "meningotype", self.samples(tmp_path)["nm_ng"], tmp_path)
        assert call(ng, "serogroup")["phenotype"] == "NG"
        g = markers(SCHEMES / "meningotype.fa")
        sample = write_sample(tmp_path, "nm_not_nm", [("nm_not_nm", g["synD"])])
        gate = typing_document(gene_dbs, "meningotype", sample, tmp_path)
        detail = call(gate, "serogroup")
        assert detail["phenotype"] == "NG"
        assert detail["notes"] == ["control gene absent"]


class TestDec:
    @staticmethod
    def genes() -> dict[str, str]:
        """The fixture markers keyed by GENE NAME — duplicate records (pic x2,
        sth x2) collapse onto the primary (first) copy, exactly like a
        screen's best-hit-per-gene fold."""
        panel: dict[str, str] = {}
        for record in iter_fasta(SCHEMES / "dec.fa"):
            gene = record.id.split("~~~")[1] if "~~~" in record.id else record.id
            panel.setdefault(gene, record.sequence)
        return panel

    @classmethod
    def samples(cls, tmp_path: Path) -> dict[str, Path]:
        g = cls.genes()
        return {
            name: write_sample(tmp_path, name, [(name, "".join(g[gene] for gene in panel))])
            for name, panel in {
                "dec_eaec": ["aggR", "pic", "uidA"],
                "dec_pic_astA": ["pic", "astA", "uidA"],
                "dec_ehec": ["stx2a", "escV", "uidA"],
                "dec_etec": ["lt", "sth", "uidA"],
                "dec_hybrid": ["stx2a", "escV", "aggR", "uidA"],
                "dec_stec_eaec": ["stx2a", "aggR", "uidA"],
                "dec_uid_only": ["uidA"],
                "dec_no_uidA": ["escV"],
                "dec_epec_atypical": ["escV", "uidA"],
            }.items()
        }

    def test_aggR_pic_calls_eaec_in_both_schemes(self, gene_dbs: Path, tmp_path: Path) -> None:
        """Given aggR + pic + the uidA control, When screened, Then both
        schemes call EAEC with high confidence (aggR satisfies even the
        strict risk-monitoring rule)."""
        document = typing_document(gene_dbs, "dec", self.samples(tmp_path)["dec_eaec"], tmp_path)
        for scheme in ("gb4789_6", "risk_monitoring"):
            assert call(document, scheme)["phenotype"] == "EAEC"
            assert call(document, scheme)["confidence"] == "high"

    def test_pic_astA_without_aggR_is_the_headline_difference(
        self, gene_dbs: Path, tmp_path: Path
    ) -> None:
        """Given pic + astA + uidA and no aggR, When screened, Then the GB
        4789.6 scheme calls EAEC (any-of aggR/pic/astA) while the
        risk-monitoring scheme falls back to non-DEC (aggR mandatory) —
        the definitional divergence between the two schemes, locked."""
        document = typing_document(
            gene_dbs, "dec", self.samples(tmp_path)["dec_pic_astA"], tmp_path
        )
        gb, risk = call(document, "gb4789_6"), call(document, "risk_monitoring")
        assert gb["phenotype"] == "EAEC"
        assert gb["confidence"] == "high"
        assert risk["phenotype"] == "non-DEC"
        assert risk["confidence"] == "low"
        assert [(c["name"], c["score"]) for c in gb["components"]] == [
            ("aggR", 0.0),
            ("pic", 1.0),
            ("astA", 1.0),
        ]

    def test_duplicate_pic_records_collapse_to_one_gene_call(
        self, gene_dbs: Path, tmp_path: Path
    ) -> None:
        """Given a sample carrying pic against a db holding two pic records,
        When screened, Then exactly one pic hit survives (blastn's
        culling_limit 1 keeps the best subject per query span, abricate
        parity) and the phenotype components list pic exactly once — the
        per-gene fold is by NAME, never by record."""
        sample = self.samples(tmp_path)["dec_pic_astA"]
        pic_hits = [
            hit
            for hit in screen_json(gene_dbs, "dec", sample)["files"][0]["hits"]
            if hit["gene"] == "pic"
        ]
        assert len(pic_hits) == 1
        assert pic_hits[0]["identity_pct"] == 100.0
        components = call(typing_document(gene_dbs, "dec", sample, tmp_path), "gb4789_6")[
            "components"
        ]
        assert [c["name"] for c in components].count("pic") == 1

    def test_stx_escV_calls_ehec_in_both_schemes(self, gene_dbs: Path, tmp_path: Path) -> None:
        """Given stx2a + escV + uidA, When screened, Then both schemes call
        EHEC (stx any-of satisfied, escV required)."""
        document = typing_document(gene_dbs, "dec", self.samples(tmp_path)["dec_ehec"], tmp_path)
        for scheme in ("gb4789_6", "risk_monitoring"):
            assert call(document, scheme)["phenotype"] == "EHEC"

    def test_lt_sth_calls_etec_in_both_schemes(self, gene_dbs: Path, tmp_path: Path) -> None:
        """Given lt + sth + uidA against a db with two sth records, When
        typed, Then both schemes call ETEC (enterotoxin any-of; the sth
        duplicate collapses to one best hit, like pic)."""
        sample = self.samples(tmp_path)["dec_etec"]
        document = typing_document(gene_dbs, "dec", sample, tmp_path)
        for scheme in ("gb4789_6", "risk_monitoring"):
            assert call(document, scheme)["phenotype"] == "ETEC"
        sth_hits = [
            hit
            for hit in screen_json(gene_dbs, "dec", sample)["files"][0]["hits"]
            if hit["gene"] == "sth"
        ]
        assert len(sth_hits) == 1

    def test_hybrid_calls_ehec_with_eaec_runner_up(self, gene_dbs: Path, tmp_path: Path) -> None:
        """Given a hybrid isolate (stx2a + escV + aggR + uidA), When
        screened, Then both schemes call EHEC with EAEC at 1.0 as the
        runner_up — the severity-order tie semantics (declaration order
        decides, the co-carried pathotype surfaces below)."""
        document = typing_document(gene_dbs, "dec", self.samples(tmp_path)["dec_hybrid"], tmp_path)
        for scheme in ("gb4789_6", "risk_monitoring"):
            detail = call(document, scheme)
            assert detail["phenotype"] == "EHEC"
            assert detail["runner_up"] == {"phenotype": "EAEC", "score": 1.0}

    def test_stx_without_escV_is_stec_with_eaec_runner_up(
        self, gene_dbs: Path, tmp_path: Path
    ) -> None:
        """Given stx2a + aggR + uidA (no escV), When screened, Then both
        schemes call STEC (stx without the LEE marker) with EAEC as the
        runner_up."""
        document = typing_document(
            gene_dbs, "dec", self.samples(tmp_path)["dec_stec_eaec"], tmp_path
        )
        for scheme in ("gb4789_6", "risk_monitoring"):
            detail = call(document, scheme)
            assert detail["phenotype"] == "STEC"
            assert detail["runner_up"]["phenotype"] == "EAEC"

    def test_uidA_only_falls_back_in_both_schemes(self, gene_dbs: Path, tmp_path: Path) -> None:
        """Given an E. coli carrying only the uidA control, When screened,
        Then both schemes output the non-DEC fallback with low
        confidence."""
        document = typing_document(
            gene_dbs, "dec", self.samples(tmp_path)["dec_uid_only"], tmp_path
        )
        for scheme in ("gb4789_6", "risk_monitoring"):
            detail = call(document, scheme)
            assert detail["phenotype"] == "non-DEC"
            assert detail["confidence"] == "low"

    def test_missing_uidA_control_gates_both_schemes(self, gene_dbs: Path, tmp_path: Path) -> None:
        """Given a marker panel without uidA, When screened, Then the
        control gate fires in both schemes: non-DEC at low confidence with
        the control note (a non-E. coli carrying DEC-like markers)."""
        document = typing_document(gene_dbs, "dec", self.samples(tmp_path)["dec_no_uidA"], tmp_path)
        for scheme in ("gb4789_6", "risk_monitoring"):
            detail = call(document, scheme)
            assert detail["phenotype"] == "non-DEC"
            assert detail["confidence"] == "low"
            assert detail["notes"] == ["control gene absent"]

    def test_escV_uidA_calls_epec_atypical_in_both_schemes(
        self, gene_dbs: Path, tmp_path: Path
    ) -> None:
        """Given escV + uidA with no bfpB and no stx, When screened, Then
        both schemes call EPEC_atypical (the fourth rule: escV required,
        bfpB and every stx subunit excluded)."""
        document = typing_document(
            gene_dbs, "dec", self.samples(tmp_path)["dec_epec_atypical"], tmp_path
        )
        for scheme in ("gb4789_6", "risk_monitoring"):
            assert call(document, scheme)["phenotype"] == "EPEC_atypical"


class TestGoldens:
    @staticmethod
    def _render(
        monkeypatch: pytest.MonkeyPatch, datadir: Path, db: str, workdir: Path, sample: str
    ) -> tuple[str, str]:
        """Screen workdir/<sample>.fa through the typed db and render both
        formats at the pinned now (chdir keeps the file path in the output
        relative, like every committed golden). The screen render is the
        pure report — designation goldens live with the typing command."""
        ensure_blast()
        (database,) = [entry for entry in discover_databases(datadir) if entry.name == db]
        assert (database.path / "typing.json").is_file()
        monkeypatch.chdir(workdir)
        report = screen_file(Path(f"{sample}.fa"), database, ScreeningParams(db=db), dbtype="nucl")
        params = ScreeningParams(db=db)
        return (
            render_json([report], params, now=PINNED_NOW),
            render_markdown([report], params, now=PINNED_NOW),
        )

    def test_doumith_4b_golden(
        self, gene_dbs: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Given the doumith fixture db and the 4b sample, When rendered as
        JSON and Markdown at the pinned now, Then each output is
        byte-identical to its committed golden."""
        sample = TestDoumith.samples(tmp_path)["lm_4b"]
        rendered, markdown = self._render(
            monkeypatch, gene_dbs, "doumith", sample.parent, sample.stem
        )
        assert rendered == (GOLDEN / "gene_doumith_4b.json").read_text(encoding="utf-8")
        assert markdown == (GOLDEN / "gene_doumith_4b.md").read_text(encoding="utf-8")

    def test_shigella_mixed_golden(
        self, gene_dbs: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Given the shigella fixture db and the two-wzx mixed sample, When
        rendered at the pinned now, Then the outputs are byte-identical to
        the committed goldens."""
        sample = TestShigella.samples(tmp_path)["sh_mixed"]
        rendered, markdown = self._render(
            monkeypatch, gene_dbs, "shigella", sample.parent, sample.stem
        )
        assert rendered == (GOLDEN / "gene_shigella_mixed.json").read_text(encoding="utf-8")
        assert markdown == (GOLDEN / "gene_shigella_mixed.md").read_text(encoding="utf-8")

    def test_dec_pic_astA_golden(
        self, gene_dbs: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Given the dec fixture db and the pic+astA headline sample, When
        rendered at the pinned now, Then the screen output is the pure
        report, byte-identical to the committed goldens — and the typing
        calls land in the gapit.typing_result/1 document (the GB 4789.6
        EAEC call beside the risk-monitoring non-DEC fallback, asserted in
        the DEC suite above)."""
        sample = TestDec.samples(tmp_path)["dec_pic_astA"]
        rendered, markdown = self._render(monkeypatch, gene_dbs, "dec", sample.parent, sample.stem)
        assert rendered == (GOLDEN / "gene_dec_pic_astA.json").read_text(encoding="utf-8")
        assert markdown == (GOLDEN / "gene_dec_pic_astA.md").read_text(encoding="utf-8")
        table = sample.parent / "dec_golden.tsv"
        screened = runner.invoke(
            app,
            [
                "screen",
                sample.name,
                "--db",
                "dec",
                "--datadir",
                str(gene_dbs),
                "-o",
                str(table),
                "--quiet",
            ],
        )
        assert screened.exit_code == 0, screened.stderr
        typed = runner.invoke(
            app,
            ["typing", str(table), "--datadir", str(gene_dbs), "--format", "json", "--quiet"],
        )
        assert typed.exit_code == 0, typed.stderr
        phenotypes = json.loads(typed.stdout)["files"][0]["phenotypes"]
        assert phenotypes["gb4789_6"]["phenotype"] == "EAEC"
        assert phenotypes["risk_monitoring"]["phenotype"] == "non-DEC"

    def test_vp_compose_golden(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Given the stage-2 serotyping fixture db (its o_group/k_group/
        serotype compose schemes) and the vp_typed sample, When rendered at
        the pinned now, Then the composed O1:K13 call is byte-identical to
        the committed goldens — the VP-style O:K composition lock."""
        datadir = tmp_path / "serotyped"
        result = runner.invoke(
            app,
            [
                "db",
                "build",
                "serotyping",
                str(TYPING / "serotyping.fa"),
                "--datadir",
                str(datadir),
                "--typing",
                str(TYPING / "typing_serotyping_v2.json"),
            ],
        )
        assert result.exit_code == 0, result.stderr
        rendered, markdown = self._render(monkeypatch, datadir, "serotyping", TYPING, "vp_typed")
        assert rendered == (GOLDEN / "gene_serotyping_vp_typed.json").read_text(encoding="utf-8")
        assert markdown == (GOLDEN / "gene_serotyping_vp_typed.md").read_text(encoding="utf-8")
