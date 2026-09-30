"""Kaptive K. pneumoniae species complex K-locus provider (cluster kind,
download-on-fetch).

Official Kaptive v3 install keyword ``kpsc_k`` (Kaptive docs, Available
databases: klebgenomics.github.io/Kaptive/db/overview.html). Kaptive v3
decentralised its databases into per-species repositories that are actively
curated on ``main``; this provider tracks the raw GenBank file directly
(URL verified 2026-09-30):

    https://raw.githubusercontent.com/klebgenomics/KpSC_surface_antigen_loci/
        main/Klebsiella_pneumoniae_Species_Complex_K.gbk

The repo also ships a ``.toml`` metadata file with upstream identity
thresholds; it is NOT fetched in v1 — a future typing.json source. The
database content is GPL-3.0, so NOTHING is bundled with gapit (MIT):
`gapit db fetch kpsc_k` downloads the gbk at fetch time and builds a
kind-cluster database (locus calls; no typing model — phenotype stays null
until a typing.json is authored). Cite Kaptive when you use the results.
"""

from gapit.providers.cluster_common import ClusterProvider

NAME = "kpsc_k"

PROVIDER = ClusterProvider(
    name=NAME,
    description="K. pneumoniae species complex K locus (Kaptive)",
    source_urls=(
        "https://raw.githubusercontent.com/klebgenomics/"
        "KpSC_surface_antigen_loci/main/"
        "Klebsiella_pneumoniae_Species_Complex_K.gbk",
    ),
    license="GPL-3.0 (database content)",
    note=(
        "klebgenomics/KpSC_surface_antigen_loci "
        "Klebsiella_pneumoniae_Species_Complex_K.gbk (Kaptive v3 keyword kpsc_k); "
        "cite Wyres et al., J Clin Microbiol 2020 (Kaptive)"
    ),
)
