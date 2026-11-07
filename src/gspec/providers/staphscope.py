"""StaphScope/Ridom spa typing provider (gene kind, download-on-fetch).

Downloads the Ridom spa reference data that StaphScope (MIT) distributes —
``sparepeats.fasta`` (865 repeat sequences, 18-30 bp each) and
``spatypes.txt`` (22842 spa types as repeat-order strings) — and rebuilds
each spa type's complete X-region sequence by concatenating its repeats in
order. Each spa type becomes one gene: the best GENE hit IS the spa type
(100% identity over the repeat array).

Screen-only: no typing model. The ``tuf`` single-copy control is a
committed wheel asset (INSDC NC_007296, public domain).

The spa reference data originates from the Ridom spaServer (academic use);
nothing is bundled with gspec (MIT) — `gspec db fetch -d sau_spa` downloads
and builds locally.
"""

import json
from collections.abc import Iterable
from pathlib import Path

from gspec.db_info_models import DbGuidance, DbReference
from gspec.fasta import iter_fasta
from gspec.providers.common import Provider
from gspec.records import Record

NAME = "sau_spa"
_RAW = (
    "https://raw.githubusercontent.com/bbeckley-hub/staphscope-typing-tool/main/"
    "staphscope/modules/spa_module/spatyper/spa_typing/main"
)
_SOURCE_FILES = ("sparepeats.fasta", "spatypes.txt")

# tuf from str_lancefield (INSDC NC_007296, public domain)
_TUF_PATH = Path(__file__).parent / "data" / "tuf.fna"


def _read_repeats(workdir: Path) -> dict[str, str]:
    repeats: dict[str, str] = {}
    for fasta in iter_fasta(workdir / "sparepeats.fasta"):
        repeats[fasta.id] = fasta.sequence.upper()
    return repeats


def transform(workdir: Path) -> Iterable[Record]:
    """Yield one Record per spa type: the concatenated repeat array.

    Records are emitted in numeric spa-type order (t001, t002, ...) — the
    order the local build was validated in. BLAST's culling_limit=1 keeps
    one hit per query window, and nested spa types (a prefix of another
    type's repeat array) tie on bitscore; the subject order therefore
    decides the survivor. Numeric order makes fetch reproducible and
    consistent with the validated local build.
    """
    repeats = _read_repeats(workdir)
    pending: list[tuple[int, str, str]] = []
    for line in (workdir / "spatypes.txt").read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("NT"):
            continue
        parts = line.split(",", 1)
        if len(parts) != 2:
            continue
        spa_name = parts[0].strip()
        if not spa_name.startswith("t") or not spa_name[1:].isdigit():
            continue
        order = [token.strip() for token in parts[1].strip().split("-")]
        if not order or not all(token.isdigit() for token in order):
            continue
        pieces: list[str] = []
        valid = True
        for token in order:
            number = int(token)
            key = f"r{number:02d}" if number < 100 else f"r{number}"
            if key not in repeats:
                key = f"r{number}"
            if key not in repeats:
                valid = False
                break
            pieces.append(repeats[key])
        if not valid or not pieces:
            continue
        pending.append((int(spa_name[1:]), spa_name, "".join(pieces)))
    for _, spa_name, sequence in sorted(pending):
        yield Record(
            db=NAME,
            gene=spa_name,
            sequence=sequence,
            accession=f"ridom:{spa_name}",
            function=("spa_type",),
            product=(
                f"S. aureus spa type {spa_name} "
                f"({len(sequence) // 24} repeat units; Ridom spaServer via StaphScope, MIT)"
            ),
            source_id=spa_name,
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
                product="S. aureus single-copy control gene tuf (EF-Tu)",
                source_id="tuf",
            )


def finalize(workdir: Path, build_dir: Path) -> None:
    """Write the per-gene floors and reads/contig screening defaults."""
    floors = {"schema": "gspec.floors/1", "default": 98.0, "genes": {"tuf": 85.0}}
    (build_dir / "floors.json").write_text(json.dumps(floors, indent=2) + "\n")
    screen = {
        "schema": "gspec.screen/1",
        "reads": {
            "min_breadth": 90.0,
            "min_identity": 95.0,
            "min_mapq": 20,
            "min_depth": 5.0,
            "min_relative_depth": 0.1,
            "relative_to": "tuf",
        },
        "contigs": {"minid": 95.0, "mincov": 85.0, "merge_fragments": True},
    }
    (build_dir / "screen.json").write_text(json.dumps(screen, indent=2) + "\n")


PROVIDER = Provider(
    name=NAME,
    description="S. aureus spa typing (Ridom spa types, screen-only)",
    vendor="Ridom spaServer via StaphScope (MIT)",
    source_urls=tuple(f"{_RAW}/{name}" for name in _SOURCE_FILES),
    dbtype="nucl",
    transform=transform,
    finalize=finalize,
    license="MIT (StaphScope); spa reference data from Ridom spaServer (academic use)",
    guidance=DbGuidance(
        references=(
            DbReference(
                citation=(
                    "Beckley B. StaphScope: a species-optimized computational"
                    " pipeline for Staphylococcus aureus genomic analysis."
                    " BMC Genomics 2026."
                ),
                doi="10.1186/s12864-026-12609-x",
                url="https://github.com/bbeckley-hub/staphscope-typing-tool",
            ),
            DbReference(
                citation=(
                    "Harmsen D, Claus H, Witte W, et al. Typing of"
                    " methicillin-resistant Staphylococcus aureus in a university"
                    " hospital setting by using novel software for spa repeat"
                    " determination and database management."
                    " J Clin Microbiol 2003;41(11):5442-8."
                ),
                doi="10.1128/JCM.41.11.5442-5448.2003",
                url="http://spaserver2.ridom.de/",
            ),
        ),
        caveats=(
            "Screen-only: the best GENE hit IS the spa type; minid>=98 for exact"
            " typing, minid>=90 for draft assemblies (flanking-region mismatch"
            " lowers coverage, not identity).",
            "Novel/unknown spa types match their closest relative — verify the"
            " repeat order against spatypes.txt before reporting.",
            "spa-negative isolates report no hits.",
            "spa reference data is from the Ridom spaServer (academic use):"
            " fetched on demand, never bundled into the gspec wheel.",
        ),
    ),
)
