"""Kaptive A. baumannii OC-locus provider (cluster kind, download-on-fetch).

Official Kaptive v3 install keyword ``ab_o`` (Kaptive docs, Available
databases: klebgenomics.github.io/Kaptive/db/overview.html) — the keyword
is ``ab_o`` while the file carries the ``OC`` (oligosaccharide cluster)
name. Kaptive v3 decentralised its databases into per-species repositories
that are actively curated on ``main``; this provider tracks the raw
GenBank file directly (URL verified 2026-09-30):

    https://raw.githubusercontent.com/johannajkenyon/
        Abaumannii_surface_polysaccharide_loci/main/Acinetobacter_baumannii_OC.gbk

The repo also ships a ``.toml`` metadata file with upstream identity
thresholds; it is NOT fetched in v1 — a future typing.json source. The
database content is GPL-3.0, so NOTHING is bundled with gapit (MIT):
`gapit db fetch ab_o` downloads the gbk at fetch time and builds a
kind-cluster database. No typing model: `gapit screen --db ab_o` reports
best-locus (OC…) calls, phenotype null. Cite Kaptive when you use results.
"""

from gapit.providers.cluster_common import ClusterProvider

NAME = "ab_o"

PROVIDER = ClusterProvider(
    name=NAME,
    description="A. baumannii OC locus — official keyword ab_o (Kaptive)",
    source_urls=(
        "https://raw.githubusercontent.com/johannajkenyon/"
        "Abaumannii_surface_polysaccharide_loci/main/Acinetobacter_baumannii_OC.gbk",
    ),
    license="GPL-3.0 (database content)",
    note=(
        "johannajkenyon/Abaumannii_surface_polysaccharide_loci "
        "Acinetobacter_baumannii_OC.gbk (Kaptive v3 keyword ab_o); "
        "cite Wyres et al., J Clin Microbiol 2020 (Kaptive)"
    ),
)
