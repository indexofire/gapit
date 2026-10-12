"""gapit.cluster/2 build-side units: gene-id resolution and the packed
PHENOTYPES cell.

``resolve_gene_ids`` is the /2 build-time transform that makes every gene
id globally unique (rules key on ids, so a shared ``/gene`` name across
loci would collide); ``colliding_genes.gbk`` is its synthetic fixture (two
records, three ``wzx`` CDS, one clashing ``/locus_tag``). The /1 tripwire
pins the kaptive_style/screening features.json bytes: /1 and untyped
builds must stay byte-identical to the committed tables. The PHENOTYPES
cell tests are pure model-level renders (no engine).
"""

from pathlib import Path

from click.testing import CliRunner

from gapit.cli import app
from gapit.cluster import BestCall, ClusterReport
from gapit.clusterbuild import resolve_gene_ids
from gapit.formats.cluster2 import (
    cluster2_tsv_file_chunk,
    cluster2_tsv_preamble,
    format_cluster2_tsv,
)
from gapit.gbfeatures import parse_genbank_features
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


_CALLED = SchemeCall(phenotype="O3/O13", score=1.0, confidence="high")
_AMBIGUOUS = SchemeCall(phenotype=None, score=0.9, confidence="ambiguous")
_FALLBACK = SchemeCall(phenotype="OUT", score=0.0, confidence="low")
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
    def test_packed_cell_in_document_order(self) -> None:
        """Given a called file with three schemes, When rendered as TSV,
        Then the PHENOTYPES cell packs name=value pairs in scheme document
        order at the /1 PHENOTYPE position."""
        report = _report({"o_group": _CALLED, "k_group": _CALLED, "serotype": _CALLED}, _BEST)

        row = cluster2_tsv_file_chunk(report, csv=False, nopath=True).rstrip("\n")

        assert row == (
            "sample.fa\tOL3_or_O13\tO3/O13\t"
            "o_group=O3/O13;k_group=O3/O13;serotype=O3/O13\t"
            "100.00\t99.50\t3\t0\t-"
        )

    def test_no_call_renders_as_dash(self) -> None:
        """Given an ambiguous (null) scheme call, When rendered, Then that
        scheme's cell entry is name=- while the called siblings keep their
        phenotypes."""
        report = _report({"o_group": _CALLED, "k_group": _AMBIGUOUS}, _BEST)

        row = cluster2_tsv_file_chunk(report, csv=False, nopath=True).rstrip("\n")

        assert "\to_group=O3/O13;k_group=-\t" in row

    def test_null_best_keeps_the_phenotypes_cell(self) -> None:
        """Given a file with no locus call, When rendered, Then the row is
        the /1 null-call dashes with the phenotypes cell still packed
        (schemes evaluate without a best locus)."""
        report = _report({"o_group": _FALLBACK, "k_group": _FALLBACK}, None)

        row = cluster2_tsv_file_chunk(report, csv=False, nopath=True).rstrip("\n")

        assert row == "sample.fa\t-\t-\to_group=OUT;k_group=OUT\t0.00\t0.00\t0\t0\t-"

    def test_header_and_csv_variant(self) -> None:
        """Given the /2 renderers, When the preamble and CSV rows build,
        Then the header carries PHENOTYPES after TYPE and the separator is
        the comma."""
        report = _report({"o_group": _CALLED}, None)

        assert cluster2_tsv_preamble(csv=False, noheader=False) == (
            "FILE\tBEST_LOCUS\tTYPE\tPHENOTYPES\tCOVERAGE\tIDENTITY\tPRESENT\tPARTIAL\tMISSING_IDS\n"
        )
        assert cluster2_tsv_preamble(csv=False, noheader=True) == ""
        assert (
            format_cluster2_tsv([report], csv=True, noheader=False, nopath=True).splitlines()[-1]
            == "sample.fa,-,-,o_group=O3/O13,0.00,0.00,0,0,-"
        )
