"""typing/2 stage 3 — the cluster-path scheme cookbook, self-checked.

Two researched cluster-database schemes as validated example documents
(``tests/data/typing/schemes/``): the Vibrio parahaemolyticus O/K panel
(van der Graaf-van Bloois et al. 2023 Kaptive databases — o_group/k_group
cluster_match rules plus the ``{o_group}:{k_group}`` compose scheme) and
the V. cholerae O1/O139 serogroup + ogawa/inaba subserotype panel (wbeT
single-SNP biology). A cluster database renders gapit.cluster/2 when its
typing document holds more than one scheme — so the vpa_ok cookbook
document installs as ONE merged database (every scheme called, the
compose rendered from its ingredients) — while single-scheme documents
(v1 or one-scheme v2) keep the byte-frozen gapit.cluster/1 surface: the
cholerae fixture installs each scheme as its own database, the Kaptive
per-scheme pattern. Samples are deterministic copies/mutations of the
fixture loci (synthetic content only).
"""

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner

from gapit.cli import app
from gapit.cluster import ClusterParams, load_features, load_typing, screen_cluster_file
from gapit.db import discover_databases
from gapit.fasta import iter_fasta
from gapit.formats.cluster import format_cluster_tsv, render_cluster_json
from gapit.formats.cluster2 import format_cluster2_tsv, render_cluster2_json
from gapit.formats.cluster_md import render_cluster2_md, render_cluster_md
from gapit.typing_cluster2 import evaluate_cluster2
from gapit.typing_engine import evaluate_typing

SCHEMES = Path(__file__).parent / "data" / "typing" / "schemes"
DATA_TYPING = Path(__file__).parent / "data" / "typing"
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
    """The vpa_ok panel as two databases (o_group and k_group schemes)."""
    datadir = tmp_path_factory.mktemp("schemes_vp")
    build_cluster_db(datadir, "vp_o", "vpa_ok.gbk", split_scheme(datadir, "vpa_ok", "o_group"))
    build_cluster_db(datadir, "vp_k", "vpa_ok.gbk", split_scheme(datadir, "vpa_ok", "k_group"))
    return datadir


@pytest.fixture(scope="module")
def vpa_ok_db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The whole three-scheme vpa_ok cookbook document as ONE merged
    database — the gapit.cluster/2 install (its gene ids are already
    unique, so build-time resolution is a no-op)."""
    datadir = tmp_path_factory.mktemp("schemes_vpa_ok")
    build_cluster_db(datadir, "vpa_ok", "vpa_ok.gbk", SCHEMES / "vpa_ok.json")
    return datadir


@pytest.fixture(scope="module")
def vc_dbs(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The cholerae panel as two databases (serogroup and subserotype
    schemes over the same wbe/wbm/wbc loci)."""
    datadir = tmp_path_factory.mktemp("schemes_vc")
    for db_name, scheme_name in (
        ("vch_serogroup", "serogroup"),
        ("vc_subserotype", "subserotype"),
    ):
        build_cluster_db(
            datadir,
            db_name,
            "cholerae.gbk",
            split_scheme(datadir, "cholerae", scheme_name),
        )
    return datadir


class TestMergedVpOkCluster2:
    def test_full_vpa_ok_document_builds_one_merged_db(self, tmp_path: Path) -> None:
        """Given the three-scheme vpa_ok document on a cluster build, When
        built, Then exit 0 — the multi-scheme install renders
        gapit.cluster/2 (the /1 one-scheme guard is gone; compose now runs
        natively on the cluster path)."""
        datadir = tmp_path / "datadir"
        datadir.mkdir()
        result = runner.invoke(
            app,
            [
                "db",
                "build",
                "vpa_ok",
                str(SCHEMES / "vpa_ok.gbk"),
                "--datadir",
                str(datadir),
                "--typing",
                str(SCHEMES / "vpa_ok.json"),
            ],
        )
        assert result.exit_code == 0, result.stderr

    def test_merged_db_screens_as_cluster2(self, vpa_ok_db: Path, tmp_path: Path) -> None:
        """Given the merged vpa_ok database and the vp_typed sample, When
        screened as JSON, Then the document is gapit.cluster/2 with the
        phenotypes keys in document order (o_group, k_group, serotype),
        serotype the composed {o}:{k} of its ingredients (the fallbacks
        flow through — the sample shares no locus with this panel)."""
        document = screen_json(vpa_ok_db, "vpa_ok", DATA_TYPING / "vp_typed.fa")
        phenotypes = document["files"][0]["phenotypes"]

        assert document["schema"] == "gapit.cluster/2"
        assert list(phenotypes) == ["o_group", "k_group", "serotype"]
        o, k = phenotypes["o_group"]["phenotype"], phenotypes["k_group"]["phenotype"]
        assert phenotypes["serotype"]["phenotype"] == f"{o}:{k}"

    def test_merged_db_best_block_carries_no_phenotype_keys(
        self, vpa_ok_db: Path, tmp_path: Path
    ) -> None:
        """Given the merged vpa_ok database and a sample carrying the
        OL3_or_O13 + KL6 loci, When screened as JSON, Then the best block
        exists, calls the locus, and carries NO phenotype keys — /2 moves
        every typing call into the phenotypes map."""
        loci = loci_of(vpa_ok_db, "vpa_ok")
        sample = write_sample(
            tmp_path,
            "vp_o3k6",
            [("ctg_o", loci["OL3_or_O13"]), ("ctg_k", loci["KL6"])],
        )
        document = screen_json(vpa_ok_db, "vpa_ok", sample)
        best = document["files"][0]["best"]

        assert best is not None
        assert best["locus"] == "OL3_or_O13"
        assert "phenotype" not in best
        assert "phenotype_detail" not in best

    def test_mixed_sample_falls_back_through_compose(self, vpa_ok_db: Path) -> None:
        """Given the vp_mixed sample (gene-path wzx/wzy fragments sharing
        no vpa_ok locus), When screened against the merged db, Then there
        is no best call — the cluster engine's no-call row, the /1 split
        twin of the unrelated-sample test — while every scheme still
        evaluates: both groups fall back low and serotype composes the
        fallback pair OUT:KUT."""
        document = screen_json(vpa_ok_db, "vpa_ok", DATA_TYPING / "vp_mixed.fa")
        file_doc = document["files"][0]

        assert file_doc["best"] is None
        assert file_doc["phenotypes"]["o_group"]["phenotype"] == "OUT"
        assert file_doc["phenotypes"]["k_group"]["phenotype"] == "KUT"
        assert file_doc["phenotypes"]["serotype"]["phenotype"] == "OUT:KUT"
        assert {call["confidence"] for call in file_doc["phenotypes"].values()} == {"low"}


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
        loci = loci_of(vc_dbs, "vch_serogroup")
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
        assert screen_json(vc_dbs, "vch_serogroup", sample)["files"][0]["best"]["phenotype"] == "O1"
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
        serogroup = screen_json(vc_dbs, "vch_serogroup", sample)["files"][0]["best"]
        assert serogroup["phenotype"] == "O1"
        assert serogroup["coverage_pct"] == 100.0
        subserotype = screen_json(vc_dbs, "vc_subserotype", sample)["files"][0]["best"]
        assert subserotype["phenotype"] == "inaba"
        assert subserotype["phenotype_detail"]["score"] == 1.0

    def test_o139_calls_o139_with_the_wbfz_note(self, vc_dbs: Path, tmp_path: Path) -> None:
        """Given the wbm locus, When screened, Then the serogroup calls
        O139 and surfaces the wbfZ junction-gene trap note."""
        sample = self.samples(vc_dbs, tmp_path)["vc_o139"]
        best = screen_json(vc_dbs, "vch_serogroup", sample)["files"][0]["best"]
        assert best["phenotype"] == "O139"
        assert "never on wbfZ alone" in best["phenotype_detail"]["notes"][0]

    def test_non_typeable_o_locus_falls_back(self, vc_dbs: Path, tmp_path: Path) -> None:
        """Given a covered non-O1/non-O139 O-antigen locus (wbc), When
        screened, Then every rule scores below cutoff (key genes missing)
        and the serogroup falls back to non-O1/non-O139."""
        sample = self.samples(vc_dbs, tmp_path)["vc_wbc"]
        best = screen_json(vc_dbs, "vch_serogroup", sample)["files"][0]["best"]
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
            typing_document.schemes[0],
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

    @pytest.mark.parametrize("sample", ["vp_o3k6", "vp_none"])
    def test_merged_vpa_ok_goldens(
        self, vpa_ok_db: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, sample: str
    ) -> None:
        """Given the merged vpa_ok database and the O3:K6 / no-coverage
        samples, When rendered in all three formats at the pinned now,
        Then each output is byte-identical to its committed cluster/2
        golden (called phenotypes vs all-fallback through compose; the
        null-best row keeps the /1 dashes with the packed PHENOTYPES
        cell)."""
        if sample == "vp_o3k6":
            loci = loci_of(vpa_ok_db, "vpa_ok")
            sample_path = write_sample(
                tmp_path,
                sample,
                [("ctg_o", loci["OL3_or_O13"]), ("ctg_k", loci["KL6"])],
            )
        else:
            sample_path = write_sample(tmp_path, sample, [("ctg", "ACGT" * 300)])
        monkeypatch.chdir(sample_path.parent)
        (database,) = [entry for entry in discover_databases(vpa_ok_db) if entry.name == "vpa_ok"]
        typing_document = load_typing(database)
        assert typing_document is not None
        params = ClusterParams(db="vpa_ok")
        report = screen_cluster_file(
            Path(f"{sample}.fa"), database, load_features(database), params
        )
        report = report.model_copy(
            update={"phenotypes": evaluate_cluster2(report, typing_document)}
        )
        assert render_cluster2_json([report], params, now=PINNED_NOW) == (
            GOLDEN / f"cluster2_{sample}.json"
        ).read_text(encoding="utf-8")
        assert render_cluster2_md([report], params, now=PINNED_NOW, document=typing_document) == (
            GOLDEN / f"cluster2_{sample}.md"
        ).read_text(encoding="utf-8")
        assert format_cluster2_tsv(
            [report], csv=False, noheader=False, nopath=False, document=typing_document
        ) == (GOLDEN / f"cluster2_{sample}.tsv").read_text(encoding="utf-8")
