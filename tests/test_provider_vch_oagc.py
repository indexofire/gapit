"""Provider metadata and committed assets for vch_oagc (V. cholerae OAGC)."""

from pathlib import Path

from gapit.db_ops import db_list_entries
from gapit.gbfeatures import parse_genbank_features
from gapit.providers.vch_oagc import PROVIDER
from gapit.typing_models import TypingDocument


def test_provider_is_fetch_on_demand_cluster_with_ncbi_source() -> None:
    assert PROVIDER.name == "vch_oagc"
    assert PROVIDER.kind == "cluster"
    assert PROVIDER.dbtype == "nucl"
    assert PROVIDER.vendor == "gapit(Hangzhou CDC)"
    assert PROVIDER.license and "INSDC" in PROVIDER.license
    assert len(PROVIDER.source_urls) == 1
    url = PROVIDER.source_urls[0]
    assert "eutils.ncbi.nlm.nih.gov" in url
    assert "LC594800" in url
    assert "LC594991" in url
    typing_path = PROVIDER.typing
    assert typing_path is not None and typing_path.is_file()
    assert PROVIDER.extra_sources


def test_committed_control_locus_and_typing_asset() -> None:
    (control,) = PROVIDER.extra_sources
    loci = parse_genbank_features(control)
    assert len(loci) == 1
    assert loci[0].features.label == "CONTROL"
    assert [gene.gene_id for gene in loci[0].features.genes] == ["ompW"]

    typing_path = PROVIDER.typing
    assert typing_path is not None
    document = TypingDocument.model_validate_json(typing_path.read_text(encoding="utf-8"))
    schemes = {scheme.name: scheme for scheme in document.schemes}
    assert set(schemes) == {"o_group", "species"}
    assert schemes["o_group"].control_gene == "ompW"
    assert len(schemes["o_group"].rules) == 190
    assert schemes["species"].rules[0].phenotype == "V. cholerae"


def test_db_list_shows_provider_metadata_and_typing(tmp_path: Path) -> None:
    entry = {item.name: item for item in db_list_entries(tmp_path)}["vch_oagc"]
    assert entry.vendor == "gapit(Hangzhou CDC)"
    assert entry.description.startswith("Vibrio cholerae O-antigen gene clusters")
    assert entry.dbtype == "nucl"
    assert entry.kind == "cluster"
    assert entry.installed is False
    assert entry.typing is True
