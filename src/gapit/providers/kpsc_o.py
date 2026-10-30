"""Kaptive K. pneumoniae species complex O-locus provider (cluster kind,
download-on-fetch).

Official Kaptive v3 install keyword ``kpsc_o`` (Kaptive docs, Available
databases: klebgenomics.github.io/Kaptive/db/overview.html). Kaptive v3
decentralised its databases into per-species repositories that are actively
curated on ``main``; this provider tracks the raw GenBank file directly
(URL verified 2026-09-30):

    https://raw.githubusercontent.com/klebgenomics/KpSC_surface_antigen_loci/
        main/Klebsiella_pneumoniae_Species_Complex_O.gbk

The repo also ships a ``.toml`` metadata file with upstream identity
thresholds; it is NOT fetched in v1 — a future typing.json source. The
database content is GPL-3.0, so NOTHING is bundled with gapit (MIT):
`gapit db fetch -d kpsc_o` downloads the gbk at fetch time and builds a
kind-cluster database. No typing model: `gapit screen --db kpsc_o` reports
best-locus (OL…) calls, phenotype null. Cite Kaptive when you use results.
"""

from gapit.db_info_models import DbGuidance, DbReference
from gapit.providers.cluster_common import ClusterProvider

NAME = "kpsc_o"

PROVIDER = ClusterProvider(
    name=NAME,
    description="K. pneumoniae species complex O locus (Kaptive)",
    vendor="Kaptive (klebgenomics)",
    source_urls=(
        "https://raw.githubusercontent.com/klebgenomics/"
        "KpSC_surface_antigen_loci/main/"
        "Klebsiella_pneumoniae_Species_Complex_O.gbk",
    ),
    license="GPL-3.0 (database content)",
    note=(
        "klebgenomics/KpSC_surface_antigen_loci "
        "Klebsiella_pneumoniae_Species_Complex_O.gbk (Kaptive v3 keyword kpsc_o); "
        "cite Wyres et al., J Clin Microbiol 2020 (Kaptive)"
    ),
    guidance=DbGuidance(
        references=(
            DbReference(
                citation=(
                    "Lam MMC, Wick RR, Judd LM, Holt KE, Wyres KL. Kaptive 2.0: updated"
                    " capsule and lipopolysaccharide locus typing for the Klebsiella"
                    " pneumoniae species complex. Microb Genom 2022;8(3):000800."
                ),
                doi="10.1099/mgen.0.000800",
            ),
        ),
        caveats=(
            "O locus types (lipopolysaccharide) are independent of K types; a complete"
            " serotype needs both.",
            "Recombination between O loci is common in K. pneumoniae, so a partial or"
            " novel locus can score ambiguously rather than falling back cleanly.",
        ),
    ),
)
