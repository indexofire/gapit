"""Kaptive K. oxytoca species complex K-locus provider (cluster kind,
download-on-fetch).

Official Kaptive v3 install keyword ``kosc_k`` (Kaptive docs, Available
databases: klebgenomics.github.io/Kaptive/db/overview.html). Kaptive v3
decentralised its databases into per-species repositories that are actively
curated on ``main``; this provider tracks the raw GenBank file directly
(URL verified 2026-09-30):

    https://raw.githubusercontent.com/klebgenomics/KoSC-surface-antigen-loci/
        main/Klebsiella_oxytoca_Species_Complex_K_locus_database.gbk

The repo also ships a ``.toml`` metadata file with upstream identity
thresholds; it is NOT fetched in v1 — a future typing.json source. The
database content is GPL-3.0, so NOTHING is bundled with gapit (MIT):
`gapit db fetch kosc_k` downloads the gbk at fetch time and builds a
kind-cluster database. No typing model: `gapit screen --db kosc_k` reports
best-locus (KL…) calls, phenotype null. Cite Kaptive when you use results.
"""

from gapit.providers.cluster_common import ClusterProvider

NAME = "kosc_k"

PROVIDER = ClusterProvider(
    name=NAME,
    description="K. oxytoca species complex K locus (Kaptive)",
    source_urls=(
        "https://raw.githubusercontent.com/klebgenomics/"
        "KoSC-surface-antigen-loci/main/"
        "Klebsiella_oxytoca_Species_Complex_K_locus_database.gbk",
    ),
    license="GPL-3.0 (database content)",
    note=(
        "klebgenomics/KoSC-surface-antigen-loci "
        "Klebsiella_oxytoca_Species_Complex_K_locus_database.gbk "
        "(Kaptive v3 keyword kosc_k); "
        "cite Wyres et al., J Clin Microbiol 2020 (Kaptive)"
    ),
)
