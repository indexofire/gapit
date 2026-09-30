"""Kaptive Klebsiella O-locus provider (cluster kind, download-on-fetch).

Same provenance as the other Kaptive providers (Wyres et al., J Clin
Microbiol 2020; klebgenomics/Kaptive): the v3 repository no longer carries
the databases, so gapit pins the archived v2.0.9 tag (verified 2026-09-30).
GPL-3.0 content — nothing is bundled; fetch downloads it. No typing model:
`gapit screen --db kaptive_o` reports best-locus (OL…) calls, phenotype null.
"""

from gapit.providers.cluster_common import ClusterProvider

NAME = "kaptive_o"

PROVIDER = ClusterProvider(
    name=NAME,
    description="Kaptive Klebsiella O antigen loci (GPL-3.0; downloaded on fetch)",
    source_urls=(
        "https://raw.githubusercontent.com/klebgenomics/Kaptive/v2.0.9/"
        "reference_database/Klebsiella_o_locus_primary_reference.gbk",
    ),
    license="GPL-3.0 (database content)",
    note=(
        "Kaptive v2.0.9 Klebsiella_o_locus_primary_reference.gbk; "
        "cite Wyres et al., J Clin Microbiol 2020 (Kaptive)"
    ),
)
