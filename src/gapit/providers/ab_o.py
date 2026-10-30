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
`gapit db fetch -d ab_o` downloads the gbk at fetch time and builds a
kind-cluster database. No typing model: `gapit screen --db ab_o` reports
best-locus (OC…) calls, phenotype null. Cite Kaptive when you use results.
"""

from gapit.db_info_models import DbGuidance, DbReference
from gapit.providers.cluster_common import ClusterProvider

NAME = "ab_o"

PROVIDER = ClusterProvider(
    name=NAME,
    description="A. baumannii OC locus — official keyword ab_o (Kaptive)",
    vendor="Kaptive (Kenyon lab)",
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
    guidance=DbGuidance(
        references=(
            DbReference(
                citation=(
                    "Cahill SM, Hall RM, Kenyon JJ. An update to the database for"
                    " Acinetobacter baumannii capsular polysaccharide locus typing extends"
                    " the extensive and diverse repertoire of genes found at and outside the"
                    " K locus. Microb Genom 2022;8(10)."
                ),
                doi="10.1099/mgen.0.000878",
                note="The A. baumannii surface-polysaccharide locus database this provider"
                " tracks (Kenyon lab), including the OC locus set.",
            ),
            DbReference(
                citation="johannajkenyon/Abaumannii_surface_polysaccharide_loci",
                url=("https://github.com/johannajkenyon/Abaumannii_surface_polysaccharide_loci"),
                note="The GenBank file actually fetched.",
            ),
        ),
        caveats=(
            "Kaptive-style locus matching is not agglutination: a locus call names the"
            " biosynthetic locus, while phase variation or mutations inside it can change"
            " the expressed antigen without changing the locus type.",
            "Novel or recombinant loci below the 95/95 floors fall back to the untypeable"
            " call rather than being forced to the nearest type.",
        ),
    ),
)
