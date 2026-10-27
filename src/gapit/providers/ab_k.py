"""Kaptive A. baumannii K-locus provider (cluster kind, download-on-fetch).

Official Kaptive v3 install keyword ``ab_k`` (Kaptive docs, Available
databases: klebgenomics.github.io/Kaptive/db/overview.html). Kaptive v3
decentralised its databases into per-species repositories that are actively
curated on ``main``; this provider tracks the raw GenBank file directly
(URL verified 2026-09-30):

    https://raw.githubusercontent.com/johannajkenyon/
        Abaumannii_surface_polysaccharide_loci/main/Acinetobacter_baumannii_K.gbk

The repo also ships a ``.toml`` metadata file with upstream identity
thresholds; it is NOT fetched in v1 — a future typing.json source. The
database content is GPL-3.0, so NOTHING is bundled with gapit (MIT):
`gapit db fetch ab_k` downloads the gbk at fetch time and builds a
kind-cluster database. No typing model: `gapit screen --db ab_k` reports
best-locus (K…) calls, phenotype null. Cite Kaptive when you use results.
"""

from gapit.db_info_models import DbGuidance, DbReference
from gapit.providers.cluster_common import ClusterProvider

NAME = "ab_k"

PROVIDER = ClusterProvider(
    name=NAME,
    description="A. baumannii K locus (Kaptive)",
    vendor="Kaptive (Kenyon lab)",
    source_urls=(
        "https://raw.githubusercontent.com/johannajkenyon/"
        "Abaumannii_surface_polysaccharide_loci/main/Acinetobacter_baumannii_K.gbk",
    ),
    license="GPL-3.0 (database content)",
    note=(
        "johannajkenyon/Abaumannii_surface_polysaccharide_loci "
        "Acinetobacter_baumannii_K.gbk (Kaptive v3 keyword ab_k); "
        "cite Wyres et al., J Clin Microbiol 2020 (Kaptive)"
    ),
    guidance=DbGuidance(
        references=(
            DbReference(
                citation=(
                    "Cahill SM, Hall RM, Kenyon JJ. An update to the database for"
                    " Acinetobacter baumannii capsular polysaccharide locus typing extends"
                    " the extensive and diverse repertoire of genes found at and outside the"
                    " K locus. Microb Genom 2022;8(10)."
                ),
                doi="10.1099/mgen.0.000878",
                note="The A. baumannii surface-polysaccharide locus database this provider"
                " tracks (Kenyon lab), including the OC locus set.",
            ),
            DbReference(
                citation="johannajkenyon/Abaumannii_surface_polysaccharide_loci",
                url=("https://github.com/johannajkenyon/Abaumannii_surface_polysaccharide_loci"),
                note="The GenBank file actually fetched.",
            ),
        ),
        caveats=(
            "Kaptive-style locus matching is not agglutination: a locus call names the"
            " biosynthetic locus, while phase variation or mutations inside it can change"
            " the expressed antigen without changing the locus type.",
            "Novel or recombinant loci below the 95/95 floors fall back to the untypeable"
            " call rather than being forced to the nearest type.",
        ),
    ),
)
