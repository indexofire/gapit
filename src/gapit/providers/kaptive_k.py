"""Kaptive Klebsiella K-locus provider (cluster kind, download-on-fetch).

Kaptive (Wyres et al., J Clin Microbiol 2020; klebgenomics/Kaptive) ships
its reference databases as GenBank locus collections. Kaptive v3
"decentralised" the databases out of the git repository (no stable raw
path), so gapit pins the ARCHIVED v2.0.9 tag of the same (transferred)
repository — stable, direct raw.githubusercontent URLs, verified 2026-09-30:

    https://raw.githubusercontent.com/klebgenomics/Kaptive/v2.0.9/
        reference_database/Klebsiella_k_locus_primary_reference.gbk

The database content is GPL-3.0, so NOTHING is bundled with gapit (MIT):
`gapit db fetch kaptive_k` downloads the gbk at fetch time and builds a
kind-cluster database (locus calls; no typing model — phenotype stays null
until a typing.json is authored). Cite Kaptive when you use the results.
"""

from gapit.providers.cluster_common import ClusterProvider

NAME = "kaptive_k"

PROVIDER = ClusterProvider(
    name=NAME,
    description="Kaptive Klebsiella K antigen loci (GPL-3.0; downloaded on fetch)",
    source_urls=(
        "https://raw.githubusercontent.com/klebgenomics/Kaptive/v2.0.9/"
        "reference_database/Klebsiella_k_locus_primary_reference.gbk",
    ),
    license="GPL-3.0 (database content)",
    note=(
        "Kaptive v2.0.9 Klebsiella_k_locus_primary_reference.gbk; "
        "cite Wyres et al., J Clin Microbiol 2020 (Kaptive)"
    ),
)
