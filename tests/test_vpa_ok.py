"""vpa_ok suite: the committed typing asset and the bundle↔build differential.

vpa_ok is a PURE built-in bundled database (the ecoli_dec/lmo_sero
pattern): no provider module, no source URLs, no `db fetch` — the wheel's
typed cluster snapshot materializes on first use and one screen renders
gapit.cluster/2. These tests pin the committed gapit.typing/2 asset's
shape and, when local copies of the real vpautils GenBank files exist
(private upstream), prove the snapshot is byte-identical to a
`gapit db build` over the merged GBKs — the regeneration contract behind
`scripts/generate_vp_typing.py` and any future snapshot refresh. The
offline bundled-path behavior (materialization, cluster/2 phenotypes)
lives in test_bundled.py; the RIMD acceptance in test_vpa_ok_acceptance.py.
"""

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

import gapit
from gapit.bundled import find_bundled
from gapit.cli import app
from gapit.fasta import iter_fasta

TYPING_ASSET = Path(gapit.__file__).resolve().parent / "data" / "dbs" / "vpa_ok" / "typing.json"
VPA_OK_LOCI = 195

runner = CliRunner()

_VPAUTILS = Path.home() / "repos/github/vpautils"
# The two real vpautils serotype references, wherever a local copy lives
# (private upstream: try the working repo, then the developer scratch dir).
REAL_OAGC = next(
    (
        candidate
        for candidate in (
            _VPAUTILS / "src/vpautils/database/serotype/OAgc.gbk",
            Path("/tmp/opencode/vp_provider_test/OAgc.gbk"),
        )
        if candidate.is_file()
    ),
    None,
)
REAL_CPSGC = next(
    (
        candidate
        for candidate in (
            _VPAUTILS / "src/vpautils/database/serotype/CPSgc.gbk",
            Path("/tmp/opencode/vp_provider_test/CPSgc.gbk"),
        )
        if candidate.is_file()
    ),
    None,
)
REAL_DATA = pytest.mark.skipif(
    REAL_OAGC is None or REAL_CPSGC is None,
    reason="local OAgc.gbk/CPSgc.gbk copies absent (private upstream)",
)


def test_committed_asset_holds_three_schemes() -> None:
    """Given the committed typing asset, When parsed, Then it is a
    gapit.typing/2 document whose serotype scheme composes the o_group
    and k_group schemes with OUT/KUT fallbacks, and every rule is a
    cluster_match over resolved (locus_tag-style) gene ids."""
    document = json.loads(TYPING_ASSET.read_text(encoding="utf-8"))
    assert document["schema"] == "gapit.typing/2"
    names = [scheme["name"] for scheme in document["schemes"]]
    assert names == ["o_group", "k_group", "serotype"]
    o_group, k_group, serotype = document["schemes"]
    assert serotype["compose"] == "{o_group}:{k_group}"
    assert serotype["fallback"] == "OUT:KUT"
    assert o_group["fallback"] == "OUT" and k_group["fallback"] == "KUT"
    for scheme in (o_group, k_group):
        assert scheme["cutoff"] == 0.9 and scheme["ambiguity_margin"] == 0.05
        for rule in scheme["rules"]:
            assert rule["model"] == "cluster_match"
            assert all("_" in gene for gene in rule["key_genes"]["genes"])


def loci_fasta(destination: Path, *loci: str) -> None:
    """Write the named loci of the bundled vpa_ok sequences as one FASTA
    record per locus — the synthetic carrier built from the wheel's own
    content, zero network."""
    bundled = find_bundled("vpa_ok")
    assert bundled is not None
    wanted = set(loci)
    with destination.open("w", encoding="utf-8") as handle:
        for record in iter_fasta(bundled.sequences_path):
            if record.id in wanted:
                handle.write(f">{record.id} sampled\n{record.sequence}\n")


@REAL_DATA
def test_bundle_is_byte_identical_to_a_db_build(tmp_path: Path) -> None:
    """Given the wheel's vpa_ok snapshot and the real vpautils GBKs it was
    produced from, When one datadir is built by `gapit db build` over the
    concatenated GBKs (--typing the committed asset) and another by
    bundled materialization (one carrier screen), Then both sequences AND
    features.json are BYTE-identical — the cluster snapshot is precisely
    what a plain build installs, gene-id resolution included."""
    assert REAL_OAGC is not None and REAL_CPSGC is not None
    merged = tmp_path / "merged.gbk"
    merged.write_bytes(REAL_OAGC.read_bytes() + REAL_CPSGC.read_bytes())

    build_dd = tmp_path / "build-dd"
    build_dd.mkdir()
    built = runner.invoke(
        app,
        [
            "db",
            "build",
            "vpa_ok",
            str(merged),
            "--typing",
            str(TYPING_ASSET),
            "--datadir",
            str(build_dd),
            "--quiet",
        ],
    )
    assert built.exit_code == 0, built.stderr
    assert json.loads(built.stdout)["records"] == VPA_OK_LOCI

    sample = tmp_path / "carrier.fna"
    loci_fasta(sample, "OL3")
    bundle_dd = tmp_path / "bundle-dd"
    bundle_dd.mkdir()
    screened = runner.invoke(
        app,
        ["screen", str(sample), "--db", "vpa_ok", "--datadir", str(bundle_dd), "--quiet"],
    )
    assert screened.exit_code == 0, screened.stderr

    for artifact in ("sequences", "features.json"):
        assert (bundle_dd / "vpa_ok" / artifact).read_bytes() == (
            build_dd / "vpa_ok" / artifact
        ).read_bytes(), artifact
