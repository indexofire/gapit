"""typing/2 stage 3 — the cluster-path scheme cookbook, self-checked.

Two researched cluster-database schemes as validated example documents
(``tests/data/typing/schemes/``): the Vibrio parahaemolyticus O/K panel
(van der Graaf-van Bloois et al. 2023 Kaptive databases — o_group/k_group
cluster_match rules plus the ``{o_group}:{k_group}`` compose scheme) and
the V. cholerae O1/O139 serogroup + ogawa/inaba subserotype panel (wbeT
single-SNP biology). A cluster database carries exactly ONE scheme
(gapit.cluster/1 has a single phenotype slot), so the runnable fixtures
install each scheme of the cookbook document as its own database — the
Kaptive O/K pattern — and the whole-document build is asserted to fail
with the typed guard. Samples are deterministic copies/mutations of the
fixture loci (synthetic content only).
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
from gapit.typing_models import single_scheme

SCHEMES = Path(__file__).parent / "data" / "typing" / "schemes"
GOLDEN = Path(__file__).parent / "golden"
PINNED_NOW = datetime(2026, 10, 1, 12, 0, 0, tzinfo=UTC)
FLIP = {"A": "C", "C": "G", "G": "T", "T": "A"}

runner = CliRunner()


def mutate(parent: str, *sites: int) -> str:
    bases = list(parent)
    for site in sites:
        bases[site] = FLIP[bases[site]]
    return "".join(bases)


def split_scheme(tmp_path: Path, document: str, scheme_name: str) -> Path:
    """One single-scheme document extracted from a cookbook file — the
    per-database split a cluster install requires."""
    document_data = json.loads((SCHEMES / f"{document}.json").read_text(encoding="utf-8"))
    scheme = next(entry for entry in document_data["schemes"] if entry["name"] == scheme_name)
    path = tmp_path / f"{document}_{scheme_name}.json"
    path.write_text(
        json.dumps({"schema": "gapit.typing/2", "schemes": [scheme]}, indent=2), encoding="utf-8"
    )
    return path


def build_cluster_db(datadir: Path, db: str, gbk: str, typing: Path) -> None:
    result = runner.invoke(
        app,
        [
            "db",
            "build",
            db,
            str(SCHEMES / gbk),
            "--datadir",
            str(datadir),
            "--typing",
            str(typing),
        ],
    )
    assert result.exit_code == 0, result.stderr


def loci_of(datadir: Path, db: str) -> dict[str, str]:
    return {record.id: record.sequence for record in iter_fasta(datadir / db / "sequences")}


def write_sample(directory: Path, name: str, contigs: list[tuple[str, str]]) -> Path:
    path = directory / f"{name}.fa"
    path.write_text(
        "".join(f">{contig}\n{sequence}\n" for contig, sequence in contigs), encoding="utf-8"
    )
    return path


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


@pytest.fixture(scope="module")
def vp_dbs(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The vp_ok panel as two databases (o_group and k_group schemes)."""
    datadir = tmp_path_factory.mktemp("schemes_vp")
    build_cluster_db(datadir, "vp_o", "vp_ok.gbk", split_scheme(datadir, "vp_ok", "o_group"))
    build_cluster_db(datadir, "vp_k", "vp_ok.gbk", split_scheme(datadir, "vp_ok", "k_group"))
    return datadir


@pytest.fixture(scope="module")
def vc_dbs(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The cholerae panel as two databases (serogroup and subserotype
    schemes over the same wbe/wbm/wbc loci)."""
    datadir = tmp_path_factory.mktemp("schemes_vc")
    for scheme_name in ("serogroup", "subserotype"):
        build_cluster_db(
            datadir,
            f"vc_{scheme_name}",
            "cholerae.gbk",
            split_scheme(datadir, "cholerae", scheme_name),
        )
    return datadir


class TestMultiSchemeClusterDocumentRejected:
    def test_full_vp_ok_document_fails_the_cluster_build(self, tmp_path: Path) -> None:
        """Given the three-scheme vp_ok document on a cluster build, When
        built, Then exit 4 with TYPING_MALFORMED naming the one-scheme
        constraint (the compose pattern runs on gene databases; cluster
        installs split the schemes)."""
        datadir = tmp_path / "datadir"
        datadir.mkdir()
        result = runner.invoke(
            app,
            [
                "db",
                "build",
                "vp_ok",
                str(SCHEMES / "vp_ok.gbk"),
                "--datadir",
                str(datadir),
                "--typing",
                str(SCHEMES / "vp_ok.json"),
            ],
        )
        assert result.exit_code == 4
        envelope = json.loads(result.stderr.splitlines()[-1])
        assert envelope["code"] == "TYPING_MALFORMED"
        assert "exactly one typing scheme" in envelope["message"]
        assert envelope["context"]["schemes"] == "o_group, k_group, serotype"


class TestVpOk:
    def test_o3_sample_calls_the_combined_label_with_note(
        self, vp_dbs: Path, tmp_path: Path
    ) -> None:
        """Given a sample carrying the O3/O13 locus, When screened against
        the o_group database, Then the call is the combined O3/O13 label
        (the loci share the same genes) with the van der Graaf citation
        note surfaced."""
        loci = loci_of(vp_dbs, "vp_o")
        sample = write_sample(tmp_path, "vp_o3", [("ctg_o", loci["OL3_or_O13"])])
        best = screen_json(vp_dbs, "vp_o", sample)["files"][0]["best"]
        assert best["phenotype"] == "O3/O13"
        assert best["phenotype_detail"]["confidence"] == "high"
        assert "O3 and O13 loci contain the same genes" in best["phenotype_detail"]["notes"][0]

    def test_o1_sample_calls_o1(self, vp_dbs: Path, tmp_path: Path) -> None:
        """Given a sample carrying the OL1 locus, When screened, Then the
        o_group call is O1 at 1.0 with the coverage/identity/key-genes
        components."""
        loci = loci_of(vp_dbs, "vp_o")
        sample = write_sample(tmp_path, "vp_o1", [("ctg_o", loci["OL1"])])
        best = screen_json(vp_dbs, "vp_o", sample)["files"][0]["best"]
        assert best["phenotype"] == "O1"
        detail = best["phenotype_detail"]
        assert detail["score"] == 1.0
        assert [component["name"] for component in detail["components"]] == [
            "coverage",
            "identity",
            "key_genes",
        ]

    def test_k6_sample_calls_k6(self, vp_dbs: Path, tmp_path: Path) -> None:
        """Given a sample carrying both antigen loci, When screened against
        the k_group database, Then the K locus calls K6 — the O:K pair the
        compose scheme documents."""
        loci = loci_of(vp_dbs, "vp_k")
        sample = write_sample(
            tmp_path,
            "vp_o3k6",
            [("ctg_o", loci_of(vp_dbs, "vp_o")["OL3_or_O13"]), ("ctg_k", loci["KL6"])],
        )
        best = screen_json(vp_dbs, "vp_k", sample)["files"][0]["best"]
        assert best["phenotype"] == "K6"

    def test_divergent_o_locus_falls_back_to_out(self, vp_dbs: Path, tmp_path: Path) -> None:
        """Given a fully covered but ~5%-divergent O locus (paired SNPs
        spaced to keep minimap2 seeds intact), When screened, Then the locus
        clears min-cluster-cov yet fails the rule's identity band and the
        scheme falls back to OUT (O untypeable)."""
        loci = loci_of(vp_dbs, "vp_o")
        degraded = mutate(
            loci["OL1"],
            *(site for base in range(20, 1400, 38) for site in (base, base + 12) if site < 1400),
        )
        sample = write_sample(tmp_path, "vp_out", [("ctg_o", degraded)])
        best = screen_json(vp_dbs, "vp_o", sample)["files"][0]["best"]
        assert best is not None
        assert best["phenotype"] == "OUT"
        assert best["phenotype_detail"]["confidence"] == "low"

    def test_unrelated_sample_makes_no_call(self, vp_dbs: Path, tmp_path: Path) -> None:
        """Given a sample with no locus coverage at all, When screened,
        Then there is no best call (the cluster engine's no-call row —
        not the fallback)."""
        sample = write_sample(tmp_path, "vp_none", [("ctg", "ACGT" * 300)])
        assert screen_json(vp_dbs, "vp_o", sample)["files"][0]["best"] is None


class TestCholerae:
    @staticmethod
    def samples(vc_dbs: Path, tmp_path: Path) -> dict[str, Path]:
        loci = loci_of(vc_dbs, "vc_serogroup")
        inaba_wbe = mutate(loci["wbe"], *[1300 + i * 5 for i in range(10)])
        return {
            "vc_ogawa": write_sample(tmp_path, "vc_ogawa", [("ctg_wbe", loci["wbe"])]),
            "vc_inaba": write_sample(tmp_path, "vc_inaba", [("ctg_wbe", inaba_wbe)]),
            "vc_o139": write_sample(tmp_path, "vc_o139", [("ctg_wbm", loci["wbm"])]),
            "vc_wbc": write_sample(tmp_path, "vc_wbc", [("ctg_wbc", loci["wbc"])]),
        }

    def test_ogawa_calls_o1_and_ogawa(self, vc_dbs: Path, tmp_path: Path) -> None:
        """Given an intact wbe locus (wbeT at 100%), When screened against
        both databases, Then the serogroup calls O1 and the subserotype
        calls ogawa through the intact wbeT determinant."""
        sample = self.samples(vc_dbs, tmp_path)["vc_ogawa"]
        assert screen_json(vc_dbs, "vc_serogroup", sample)["files"][0]["best"]["phenotype"] == "O1"
        detail = screen_json(vc_dbs, "vc_subserotype", sample)["files"][0]["best"][
            "phenotype_detail"
        ]
        assert "Ogawa/Inaba determinant is single mutations in the wbeT" in detail["notes"][0]

    def test_inaba_calls_o1_and_inaba(self, vc_dbs: Path, tmp_path: Path) -> None:
        """Given the wbe locus with a mutated wbeT (98% identity, below the
        99.9 allele floor but above the locus band), When screened, Then
        the serogroup still calls O1 while the subserotype calls inaba
        through the negative-wbeT rule."""
        sample = self.samples(vc_dbs, tmp_path)["vc_inaba"]
        serogroup = screen_json(vc_dbs, "vc_serogroup", sample)["files"][0]["best"]
        assert serogroup["phenotype"] == "O1"
        assert serogroup["coverage_pct"] == 100.0
        subserotype = screen_json(vc_dbs, "vc_subserotype", sample)["files"][0]["best"]
        assert subserotype["phenotype"] == "inaba"
        assert subserotype["phenotype_detail"]["score"] == 1.0

    def test_o139_calls_o139_with_the_wbfz_note(self, vc_dbs: Path, tmp_path: Path) -> None:
        """Given the wbm locus, When screened, Then the serogroup calls
        O139 and surfaces the wbfZ junction-gene trap note."""
        sample = self.samples(vc_dbs, tmp_path)["vc_o139"]
        best = screen_json(vc_dbs, "vc_serogroup", sample)["files"][0]["best"]
        assert best["phenotype"] == "O139"
        assert "never on wbfZ alone" in best["phenotype_detail"]["notes"][0]

    def test_non_typeable_o_locus_falls_back(self, vc_dbs: Path, tmp_path: Path) -> None:
        """Given a covered non-O1/non-O139 O-antigen locus (wbc), When
        screened, Then every rule scores below cutoff (key genes missing)
        and the serogroup falls back to non-O1/non-O139."""
        sample = self.samples(vc_dbs, tmp_path)["vc_wbc"]
        best = screen_json(vc_dbs, "vc_serogroup", sample)["files"][0]["best"]
        assert best["locus"] == "wbc"
        assert best["phenotype"] == "non-O1/non-O139"
        assert best["phenotype_detail"]["score"] < 0.9

    def test_o139_sample_leaves_the_subserotype_undetermined(
        self, vc_dbs: Path, tmp_path: Path
    ) -> None:
        """Given the wbm locus against the subserotype database, When
        screened, Then no wbe gene is present and the scheme falls back to
        undetermined (subserotyping applies to O1 only)."""
        sample = self.samples(vc_dbs, tmp_path)["vc_o139"]
        best = screen_json(vc_dbs, "vc_subserotype", sample)["files"][0]["best"]
        assert best["phenotype"] == "undetermined"


class TestGoldens:
    def test_vc_inaba_goldens(
        self, vc_dbs: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Given the vc_subserotype database and the inaba sample, When
        rendered in all three formats at the pinned now, Then each output is
        byte-identical to its committed golden."""
        sample = TestCholerae.samples(vc_dbs, tmp_path)["vc_inaba"]
        monkeypatch.chdir(sample.parent)
        (database,) = [
            entry for entry in discover_databases(vc_dbs) if entry.name == "vc_subserotype"
        ]
        typing_document = load_typing(database)
        assert typing_document is not None
        params = ClusterParams(db="vc_subserotype")
        report = evaluate_typing(
            screen_cluster_file(Path("vc_inaba.fa"), database, load_features(database), params),
            single_scheme(typing_document),
        )
        assert render_cluster_json([report], params, now=PINNED_NOW, typed=True) == (
            GOLDEN / "cluster_vc_inaba.json"
        ).read_text(encoding="utf-8")
        assert render_cluster_md([report], params, now=PINNED_NOW, typed=True) == (
            GOLDEN / "cluster_vc_inaba.md"
        ).read_text(encoding="utf-8")
        assert format_cluster_tsv(
            [report], csv=False, noheader=False, nopath=False, typed=True
        ) == (GOLDEN / "cluster_vc_inaba.tsv").read_text(encoding="utf-8")
