"""gapit.cluster/2 build-side units: gene-id resolution and the per-scheme
winner-evidence columns.

``resolve_gene_ids`` is the /2 build-time transform that makes every gene
id globally unique (rules key on ids, so a shared ``/gene`` name across
loci would collide); ``colliding_genes.gbk`` is its synthetic fixture (two
records, three ``wzx`` CDS, one clashing ``/locus_tag``). The /1 tripwire
pins the kaptive_style/screening features.json bytes: /1 and untyped
builds must stay byte-identical to the committed tables. The column tests
are pure model-level renders over the vp_ok cookbook document (no engine).
"""

from pathlib import Path

from click.testing import CliRunner

from gapit.cli import app
from gapit.cluster import BestCall, ClusterReport
from gapit.clusterbuild import resolve_gene_ids
from gapit.formats.cluster2 import (
    cluster2_columns,
    cluster2_tsv_file_chunk,
    cluster2_tsv_preamble,
    format_cluster2_tsv,
)
from gapit.gbfeatures import parse_genbank_features
from gapit.typing_models import TypingDocument, read_typing_document
from gapit.typing_results import SchemeCall

DATA = Path(__file__).parent / "data" / "cluster"
GOLDEN = Path(__file__).parent / "golden"

runner = CliRunner()


class TestResolveGeneIds:
    def test_colliding_genes_resolve_to_unique_ids(self) -> None:
        """Given the colliding fixture (wzx on both records, three CDS),
        When resolved, Then unique-count ids stay verbatim (manC, gnd), a
        colliding gene takes its free /locus_tag (COL1_wzx), a colliding
        gene whose locus_tag is already claimed synthesizes
        locus:id@start, and every id is globally unique."""
        loci = parse_genbank_features(DATA / "colliding_genes.gbk")

        resolved = resolve_gene_ids(loci)

        assert [[gene.gene_id for gene in locus.features.genes] for locus in resolved] == [
            ["COL1_wzx", "manC"],
            ["COL2:wzx@101", "COL2:wzx@21", "gnd"],
        ]
        ids = [gene.gene_id for locus in resolved for gene in locus.features.genes]
        assert len(ids) == len(set(ids))

    def test_declaration_order_survives_coordinate_resolution(self) -> None:
        """Given COL2's out-of-coordinate-order CDS declarations, When
        resolved, Then ids are assigned by start coordinate (the @21 gene
        claims nothing the @101 gene needed) while the genes tuple keeps
        its declaration order — features render in file order."""
        (locus,) = resolve_gene_ids(parse_genbank_features(DATA / "colliding_genes.gbk"))[1:]

        assert [(gene.gene_id, gene.start) for gene in locus.features.genes] == [
            ("COL2:wzx@101", 101),
            ("COL2:wzx@21", 21),
            ("gnd", 181),
        ]

    def test_unique_only_loci_pass_through_verbatim(self) -> None:
        """Given a fixture with no id collisions (the vp_ok panel), When
        resolved, Then every gene keeps its primary id."""
        loci = parse_genbank_features(
            Path(__file__).parent / "data" / "typing" / "schemes" / "vp_ok.gbk"
        )

        resolved = resolve_gene_ids(loci)

        assert [gene.gene_id for locus in resolved for gene in locus.features.genes] == [
            gene.gene_id for locus in loci for gene in locus.features.genes
        ]


class TestV1FeatureBytesTripwire:
    def test_kaptive_style_features_stay_byte_identical(self, tmp_path: Path) -> None:
        """Given the kaptive-style /1 fixture, When built untyped, Then its
        features.json is byte-identical to the committed table — the /1
        tripwire (the additive GeneFeature.locus_tag never serializes)."""
        result = runner.invoke(
            app,
            ["db", "build", "tw", str(DATA / "kaptive_style.gbk"), "--datadir", str(tmp_path)],
        )
        assert result.exit_code == 0, result.stderr
        assert (tmp_path / "tw" / "features.json").read_text(encoding="utf-8") == (
            GOLDEN / "features_kaptive_style.json"
        ).read_text(encoding="utf-8")
        assert "locus_tag" not in (tmp_path / "tw" / "features.json").read_text(encoding="utf-8")

    def test_screening_features_stay_byte_identical(self, tmp_path: Path) -> None:
        """Given the screening /1 fixture (single-scheme typed builds ride
        the same path), When built, Then its features.json is
        byte-identical to the committed table."""
        result = runner.invoke(
            app,
            ["db", "build", "tw", str(DATA / "screening.gbk"), "--datadir", str(tmp_path)],
        )
        assert result.exit_code == 0, result.stderr
        assert (tmp_path / "tw" / "features.json").read_text(encoding="utf-8") == (
            GOLDEN / "features_screening.json"
        ).read_text(encoding="utf-8")


def _report(phenotypes: dict[str, SchemeCall], best: BestCall | None) -> ClusterReport:
    return ClusterReport(file="sample.fa", best=best, loci=(), phenotypes=phenotypes)


_CALLED = SchemeCall(
    phenotype="O3/O13",
    score=1.0,
    confidence="high",
    winner_locus="OL3_or_O13",
    winner_coverage_pct=100.0,
    winner_identity_pct=99.5,
    winner_present=3,
    winner_partial=0,
    winner_missing=(),
)
_K_CALLED = SchemeCall(
    phenotype="K6",
    score=1.0,
    confidence="high",
    winner_locus="KL6",
    winner_coverage_pct=100.0,
    winner_identity_pct=99.5,
    winner_present=3,
    winner_partial=0,
    winner_missing=(),
)
_SEROTYPE = SchemeCall(phenotype="O3/O13:K6", score=1.0, confidence="high")
_AMBIGUOUS = SchemeCall(phenotype=None, score=0.9, confidence="ambiguous")
_FALLBACK = SchemeCall(phenotype="OUT", score=0.0, confidence="low")
_K_FALLBACK = SchemeCall(phenotype="KUT", score=0.0, confidence="low")
_DOCUMENT = read_typing_document(
    Path(__file__).parent / "data" / "typing" / "schemes" / "vp_ok.json"
)
_BEST = BestCall(
    locus="OL3_or_O13",
    label="OL3",
    type="O3/O13",
    coverage_pct=100.0,
    identity_pct=99.5,
    genes_present=3,
    genes_partial=0,
    genes_absent=0,
)


class TestPhenotypesCell:
    def test_columns_put_compose_first_then_scheme_blocks(self) -> None:
        """Given the vp_ok cookbook document (o_group, k_group, serotype
        compose), When the /2 columns derive, Then SEROTYPE leads and each
        non-compose scheme carries its call + winner-evidence block (the
        o_group prefix shortens to O)."""
        assert cluster2_columns(_DOCUMENT) == (
            "FILE",
            "SEROTYPE",
            "O_GROUP",
            "O_LOCUS",
            "O_COVERAGE",
            "O_IDENTITY",
            "O_PRESENT",
            "O_PARTIAL",
            "O_MISSING",
            "K_GROUP",
            "K_LOCUS",
            "K_COVERAGE",
            "K_IDENTITY",
            "K_PRESENT",
            "K_PARTIAL",
            "K_MISSING",
        )

    def test_row_renders_compose_then_calls_with_evidence(self) -> None:
        """Given a called file with both schemes carrying winner evidence,
        When rendered as TSV, Then the serotype value leads and each
        scheme's block shows its winning locus's own numbers."""
        report = _report({"o_group": _CALLED, "k_group": _K_CALLED, "serotype": _SEROTYPE}, _BEST)

        row = cluster2_tsv_file_chunk(report, csv=False, nopath=True, document=_DOCUMENT).rstrip(
            "\n"
        )

        assert row == (
            "sample.fa\tO3/O13:K6\t"
            "O3/O13\tOL3_or_O13\t100.00\t99.50\t3\t0\t-\t"
            "K6\tKL6\t100.00\t99.50\t3\t0\t-"
        )

    def test_no_call_renders_as_dash(self) -> None:
        """Given an ambiguous (null) scheme call without winner evidence,
        When rendered, Then that scheme's call cell is ``-`` and its
        evidence cells follow as dashes while the called sibling keeps its
        block."""
        report = _report(
            {
                "o_group": _CALLED,
                "k_group": _AMBIGUOUS,
                "serotype": SchemeCall(phenotype=None, score=0.9, confidence="ambiguous"),
            },
            _BEST,
        )

        row = cluster2_tsv_file_chunk(report, csv=False, nopath=True, document=_DOCUMENT).rstrip(
            "\n"
        )

        assert "\tK6\t" not in row
        assert row == (
            "sample.fa\t-\tO3/O13\tOL3_or_O13\t100.00\t99.50\t3\t0\t-\t-\t-\t-\t-\t-\t-\t-"
        )

    def test_fallback_calls_keep_dashed_evidence(self) -> None:
        """Given a file with no locus call and both schemes at fallback,
        When rendered, Then each call column carries its fallback label
        with dashed evidence (schemes evaluate without a best locus)."""
        report = _report(
            {
                "o_group": _FALLBACK,
                "k_group": _K_FALLBACK,
                "serotype": SchemeCall(phenotype="OUT:KUT", score=0.0, confidence="low"),
            },
            None,
        )

        row = cluster2_tsv_file_chunk(report, csv=False, nopath=True, document=_DOCUMENT).rstrip(
            "\n"
        )

        assert row == ("sample.fa\tOUT:KUT\tOUT\t-\t-\t-\t-\t-\t-\tKUT\t-\t-\t-\t-\t-\t-")

    def test_header_and_csv_variant(self) -> None:
        """Given the /2 renderers, When the preamble and CSV rows build,
        Then the header derives from the document and the separator is the
        comma (a KUT sibling shows the mixed evidence/dash shape)."""
        report = _report(
            {
                "o_group": _CALLED,
                "k_group": _K_FALLBACK,
                "serotype": SchemeCall(phenotype="O3/O13:KUT", score=0.9, confidence="high"),
            },
            None,
        )

        assert cluster2_tsv_preamble(csv=False, noheader=False, document=_DOCUMENT) == (
            "\t".join(cluster2_columns(_DOCUMENT)) + "\n"
        )
        assert cluster2_tsv_preamble(csv=False, noheader=True, document=_DOCUMENT) == ""
        assert (
            format_cluster2_tsv(
                [report], csv=True, noheader=False, nopath=True, document=_DOCUMENT
            ).splitlines()[-1]
            == "sample.fa,O3/O13:KUT,O3/O13,OL3_or_O13,100.00,99.50,3,0,-,KUT,-,-,-,-,-,-"
        )


def test_colliding_evidence_prefixes_fall_back_to_full_tokens() -> None:
    """Given schemes o_group and o (both preferring O), When /2 columns are
    derived, Then they use O_GROUP/O prefixes so no column name repeats."""
    rule = {
        "model": "weighted_genes",
        "phenotype": "X",
        "weights": {"g": 1.0},
        "identity_floor": 90.0,
    }
    document = TypingDocument.model_validate(
        {
            "schemes": [
                {
                    "name": "o_group",
                    "rules": [rule],
                    "cutoff": 0.5,
                    "ambiguity_margin": 0.1,
                    "fallback": "OUT",
                },
                {
                    "name": "o",
                    "rules": [rule],
                    "cutoff": 0.5,
                    "ambiguity_margin": 0.1,
                    "fallback": "OUT",
                },
            ]
        }
    )

    columns = cluster2_columns(document)

    assert len(columns) == len(set(columns))
    assert "O_GROUP_LOCUS" in columns
    assert "O_LOCUS" in columns
