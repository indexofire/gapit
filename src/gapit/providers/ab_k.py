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
                    "Wyres KL, Wick RR, Gorrie C, Jenney A, Follador R, Thomson NR, Holt KE."
                    " Identification of Klebsiella capsule synthesis loci from whole genome"
                    " data. Microb Genom 2016;2(12):e000102."
                ),
                doi="10.1099/mgen.0.000102",
                note="The Kaptive method this provider's locus matching follows.",
            ),
        ),
        caveats=(
            "K locus typing describes the capsule biosynthesis locus, not the expressed"
            " capsule: phase variation and mutations inside the locus can change the"
            " serotype without changing the locus type.",
            "Kaptive-style scoring is not agglutination; novel or recombinant loci below"
            " the 95/95 floors fall back to the untypeable call.",
        ),
    ),
)
