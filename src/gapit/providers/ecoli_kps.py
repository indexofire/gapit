"""Kaptive E. coli group 2+3 capsular polysaccharide provider (cluster
kind, download-on-fetch).

Official Kaptive v3 install keyword ``ecoli_kps`` (Kaptive docs, Available
databases: klebgenomics.github.io/Kaptive/db/overview.html). Kaptive v3
decentralised its databases into per-species repositories that are actively
curated on ``main``; this provider tracks the raw GenBank file directly
(URL verified 2026-09-30):

    https://raw.githubusercontent.com/rgladstone/EC-K-typing/main/
        EC-K-typing_group2and3.gbk

The repo also ships a ``.toml`` metadata file with upstream identity
thresholds; it is NOT fetched in v1 — a future typing.json source. The
database content is GPL-3.0, so NOTHING is bundled with gapit (MIT):
`gapit db fetch ecoli_kps` downloads the gbk at fetch time and builds a
kind-cluster database. No typing model: `gapit screen --db ecoli_kps`
reports best-locus (group 2+3 K…) calls, phenotype null. Cite Kaptive when
you use the results.
"""

from gapit.providers.cluster_common import ClusterProvider

NAME = "ecoli_kps"

PROVIDER = ClusterProvider(
    name=NAME,
    description="E. coli group 2+3 capsular polysaccharide loci (Kaptive)",
    source_urls=(
        "https://raw.githubusercontent.com/rgladstone/EC-K-typing/main/EC-K-typing_group2and3.gbk",
    ),
    license="GPL-3.0 (database content)",
    note=(
        "rgladstone/EC-K-typing EC-K-typing_group2and3.gbk "
        "(Kaptive v3 keyword ecoli_kps); "
        "cite Wyres et al., J Clin Microbiol 2020 (Kaptive)"
    ),
)
