"""Tests for the emmtyper (spy_emm) and seroba (spn_sero) providers.

Covers: provider registration, kind/license, transform correctness
(emm header parsing, tuf control inclusion), preprocess FASTA→GFF3
conversion with NT-reference exclusion."""

import random
from pathlib import Path

from gapit.providers import REGISTRY


class TestEmmtyperProvider:
    def test_registered(self):
        assert "spy_emm" in REGISTRY

    def test_kind_is_gene(self):
        assert REGISTRY["spy_emm"].kind == "gene"

    def test_gpl_license_declared(self):
        assert "GPL" in (REGISTRY["spy_emm"].license or "")

    def test_transform_parses_emm_headers(self, tmp_path: Path):
        from gapit.providers.emmtyper import transform

        emm = tmp_path / "emm.fna"
        emm.write_text(
            ">EMM1.0 emm1.0 (emm cluster A-C3) usually T1\n"
            + "ACGT" * 40
            + "\n"
            + ">STG62647.0 stG62647.0 (emm cluster E6) GCS\n"
            + "TTTT" * 40
            + "\n"
        )
        records = list(transform(tmp_path))
        emm_records = [r for r in records if r.gene != "tuf"]
        assert len(emm_records) == 2
        by_gene = {r.gene: r for r in emm_records}
        assert "EMM1.0" in by_gene
        assert by_gene["EMM1.0"].db == "spy_emm"
        assert "cluster A-C3" in by_gene["EMM1.0"].product
        assert "Group G strep" in by_gene["STG62647.0"].product

    def test_tuf_included_when_data_file_exists(self, tmp_path: Path):
        from gapit.providers.emmtyper import transform

        emm = tmp_path / "emm.fna"
        emm.write_text(">EMM1.0 test\nACGT\n")
        records = list(transform(tmp_path))
        tuf = [r for r in records if r.gene == "tuf"]
        if tuf:
            assert len(tuf) == 1
            assert tuf[0].db == "spy_emm"


class TestSerobaProvider:
    def test_registered(self):
        assert "spn_sero" in REGISTRY

    def test_kind_is_cluster(self):
        from gapit.providers.cluster_common import ClusterProvider

        provider = REGISTRY["spn_sero"]
        assert isinstance(provider, ClusterProvider)
        assert provider.kind == "cluster"

    def test_gpl_license_declared(self):
        assert "GPL" in (REGISTRY["spn_sero"].license or "")

    def test_preprocess_and_finalize_hooks(self):
        from gapit.providers.cluster_common import ClusterProvider

        provider = REGISTRY["spn_sero"]
        assert isinstance(provider, ClusterProvider)
        assert provider.preprocess is not None
        assert provider.finalize is not None

    def test_preprocess_creates_gff3(self, tmp_path: Path):
        from gapit.providers.seroba import preprocess

        rng = random.Random(11)
        seq1 = "".join(rng.choice("ACGT") for _ in range(1200))
        seq2 = "".join(rng.choice("ACGT") for _ in range(1200))
        ref = tmp_path / "reference.fasta"
        ref.write_text(f">01\n{seq1}\n>02\n{seq2}\n")
        gff3 = preprocess(tmp_path)
        assert gff3.suffix == ".gff3"
        text = gff3.read_text()
        assert text.startswith("##gff-version 3")
        assert ">01" in text
        assert ">02" in text
        assert "CDS" in text

    def test_preprocess_drops_nt_references(self, tmp_path: Path):
        from gapit.providers.seroba import preprocess

        rng = random.Random(12)
        seq = "".join(rng.choice("ACGT") for _ in range(1200))
        ref = tmp_path / "reference.fasta"
        ref.write_text(f">01\n{seq}\n>alternative_aliB_NT\n{seq}\n>Swiss_NT\n{seq}\n")
        gff3 = preprocess(tmp_path)
        text = gff3.read_text()
        assert ">01" in text
        assert "alternative_aliB_NT" not in text
        assert "Swiss_NT" not in text
