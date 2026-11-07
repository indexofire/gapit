"""emmtyper CDC emm database provider (gene kind, download-on-fetch).

Upstream is the CDC emm allele database distributed via emmtyper (MDU-PHL,
GPLv3). The provider downloads ``emm.fna`` (2113 references, 5' hypervariable
region, 156-183 bp each) from the emmtyper GitHub repository and builds a
gene database with dual typing schemes: emm_type (356 rules, identity_floor
98% to filter emm-like cross-reactions from mrp/enn at ~97%) and emm_cluster
(50 rules, 95% floor). The tuf single-copy control gene is a committed
wheel asset (from str_lancefield's INSDC-sourced tuf).

GPL-3.0: nothing is bundled with gspec (MIT); `gspec db fetch spy_emm`
downloads and builds locally. Cite CDC/emmtyper when using results.
"""

import json
import re
from collections.abc import Iterable
from pathlib import Path

from gapit.db_info_models import DbGuidance, DbReference
from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record

NAME = "spy_emm"
_RAW = "https://raw.githubusercontent.com/MDU-PHL/emmtyper/master/emmtyper/db"
_SOURCE_FILES = ("emm.fna",)

# tuf from str_lancefield (INSDC NC_007296, public domain)
_TUF_PATH = Path(__file__).parent / "data" / "tuf.fna"


def transform(workdir: Path) -> Iterable[Record]:
    """Yield emm allele Records from ``workdir/emm.fna`` plus the tuf control."""
    for fasta in iter_fasta(workdir / "emm.fna"):
        allele = fasta.id  # e.g., EMM1.0, STG62647.0
        cluster = ""
        match = re.search(r"cluster ([A-Z]\d*(?:-[A-Z]\d*)*)", fasta.description)
        if match:
            cluster = match.group(1)
        product = f"emm allele {allele}"
        if cluster:
            product += f" (cluster {cluster})"
        if allele.startswith("STG"):
            product += "; Group G strep emm-like (stG)"
        elif allele.startswith("STC"):
            product += "; Group C strep emm-like (stC)"
        elif allele.startswith("STL"):
            product += "; Group L strep emm-like (stL)"
        else:
            product += "; S. pyogenes emm"
        yield Record(
            db=NAME,
            gene=allele,
            sequence=fasta.sequence,
            accession=f"emmtyper:{allele}",
            function=("serotype",),
            product=product,
            source_id=fasta.id,
        )
    # tuf control (committed wheel asset)
    if _TUF_PATH.is_file():
        for fasta in iter_fasta(_TUF_PATH):
            yield Record(
                db=NAME,
                gene="tuf",
                sequence=fasta.sequence,
                accession="str_lancefield:tuf",
                function=("control",),
                product="S. pyogenes single-copy control gene tuf (EF-Tu)",
                source_id=fasta.id,
            )


def finalize(workdir: Path, build_dir: Path) -> None:
    """Write typing, screen, and floors sidecars from committed assets."""
    typing_path = Path(__file__).parent / "data" / "spy_emm_typing.json"
    if typing_path.is_file():
        (build_dir / "typing.json").write_text(
            typing_path.read_text(encoding="utf-8"), encoding="utf-8"
        )
    screen = {
        "schema": "gapit.screen/1",
        "reads": {
            "min_breadth": 90.0,
            "min_identity": 90.0,
            "min_mapq": 20,
            "min_depth": 5.0,
            "min_relative_depth": 0.1,
            "relative_to": "tuf",
        },
        "contigs": {"minid": 95.0, "mincov": 85.0, "merge_fragments": True},
    }
    (build_dir / "screen.json").write_text(json.dumps(screen, indent=2) + "\n")
    floors = {"schema": "gapit.floors/1", "default": 95.0, "genes": {"tuf": 85.0}}
    (build_dir / "floors.json").write_text(json.dumps(floors, indent=2) + "\n")


PROVIDER = Provider(
    name=NAME,
    description="S. pyogenes emm typing (CDC alleles, type + cluster)",
    vendor="CDC/emmtyper (MDU-PHL)",
    source_urls=tuple(f"{_RAW}/{name}" for name in _SOURCE_FILES),
    dbtype="nucl",
    transform=transform,
    finalize=finalize,
    license="GPL-3.0 (emmtyper/CDC emm database)",
    guidance=DbGuidance(
        references=(
            DbReference(
                citation=(
                    "McMillan DJ, Dreze PA, Vu T, et al. Updated model of group A"
                    " Streptococcus M proteins based on a comprehensive worldwide study."
                    " Clin Microbiol Infect 2013;19(5):E222-9."
                ),
                doi="10.1111/1469-0691.12149",
            ),
            DbReference(
                citation=(
                    "Sanderson-Smith M, De Oliveira DM, Guglielmini J, et al. A systematic"
                    " and functional classification of S. pyogenes serves as a new tool for"
                    " molecular typing and vaccine development. J Infect Dis 2014;209(8):1325-38."
                ),
                doi="10.1093/infdis/jit643",
            ),
        ),
        caveats=(
            "Screen GENES column shows the exact subtype allele (e.g., EMM28.0);"
            " typing adds type/cluster labels",
            "emm_type identity_floor 98% filters emm-like cross-reactions (mrp/enn"
            " ~97% identity); subtypes with <98% identity to the reference may be missed",
            "GAS strains carry 1-3 emm-family genes; only the true emm allele reaches"
            " 100% identity",
            "GPL-3.0: fetched on demand, never bundled into the gspec wheel",
        ),
    ),
)
