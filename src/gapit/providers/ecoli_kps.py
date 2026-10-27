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

from gapit.db_info_models import DbGuidance, DbReference
from gapit.providers.cluster_common import ClusterProvider

NAME = "ecoli_kps"

PROVIDER = ClusterProvider(
    name=NAME,
    description="E. coli group 2+3 capsular polysaccharide loci (Kaptive)",
    vendor="Kaptive (Gladstone lab)",
    source_urls=(
        "https://raw.githubusercontent.com/rgladstone/EC-K-typing/main/EC-K-typing_group2and3.gbk",
    ),
    license="GPL-3.0 (database content)",
    note=(
        "rgladstone/EC-K-typing EC-K-typing_group2and3.gbk "
        "(Kaptive v3 keyword ecoli_kps); "
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
                note="Kaptive method; this provider tracks the EC-K-typing group 2+3 locus set.",
            ),
        ),
        caveats=(
            "Only E. coli group 2 and group 3 capsules are covered; group 1 and group 4"
            " capsules have no locus here and fall back to the untypeable call.",
            "The K locus is not the O antigen and is not the pathotype: K typing answers"
            " one capsule question only.",
        ),
    ),
)
