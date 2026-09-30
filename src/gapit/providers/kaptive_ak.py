"""Kaptive Acinetobacter baumannii K-locus provider (cluster kind,
download-on-fetch).

Same provenance as the other Kaptive providers (Wyres et al., J Clin
Microbiol 2020; klebgenomics/Kaptive): the v3 repository no longer carries
the databases, so gapit pins the archived v2.0.9 tag (verified 2026-09-30).
GPL-3.0 content — nothing is bundled; fetch downloads it. No typing model:
`gapit screen --db kaptive_ak` reports best-locus (K…) calls, phenotype null.
"""

from gapit.providers.cluster_common import ClusterProvider

NAME = "kaptive_ak"

PROVIDER = ClusterProvider(
    name=NAME,
    description="Kaptive A. baumannii K antigen loci (GPL-3.0; downloaded on fetch)",
    source_urls=(
        "https://raw.githubusercontent.com/klebgenomics/Kaptive/v2.0.9/"
        "reference_database/Acinetobacter_baumannii_k_locus_primary_reference.gbk",
    ),
    license="GPL-3.0 (database content)",
    note=(
        "Kaptive v2.0.9 Acinetobacter_baumannii_k_locus_primary_reference.gbk; "
        "cite Wyres et al., J Clin Microbiol 2020 (Kaptive)"
    ),
)
