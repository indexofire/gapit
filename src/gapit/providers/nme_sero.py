"""nme_sero provider: Neisseria meningitidis serogroup/ctrA/sodC panel.

The non-MLST, gene-presence part of MDU-PHL/meningotype (GPL-3.0): the
serogroup A/B/C/W/X/Y amplicons, the ctrA capsule marker, and a single-copy
sodC depth control. MLST and allele-level finetyping (PorA/FetA/PorB/BAST)
are deliberately absent: MLST/cc belongs to a gMLST pipeline, and allele
calling is outside gapit's gene-presence contract.

The gapit-native nucleotide database carries the eight records, so
``gapit screen -d nme_sero`` and ``gapit typing`` work from FASTQ or
contigs. sodC is the reads-mode relative-depth control (``screen.json`` sets
``relative_to: sodC``); its reference is NC_013016.1:1293529-1294089
(N. meningitidis alpha14), embedded here because the upstream meningotype
repository does not ship a sodC FASTA.
"""

from __future__ import annotations

import json
from collections.abc import Iterator
from pathlib import Path

from gapit.db_info_models import DbGuidance, DbReference
from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record

NAME = "nme_sero"
RAW = "https://raw.githubusercontent.com/MDU-PHL/meningotype/master/meningotype/db"

_SODC = (
    "AATGAATATGAAAACCTTATTAGCACTAGCGGTTAGTGCAGTATGTTCAGTTAGTGTTGCGCAAGCACAC"
    "GAGCATAATACGATACCTAAAGGTGCTTCTATTGAAGTGAAAGTGCAACAACTTGATCCAGTAAACGGTA"
    "ACAAAGATGTGGGTACAGTGACTATTACTGAATCTAACTATGGTCTTGTGTTTACCCCTGATTTACAAGG"
    "ATTAAGCGAAGGCTTACATGGTTTCCACATCCATGAAAACCCAAGCTGTGAGCCAAAAGAAAAAGAAGGT"
    "AAATTGACAGCTGGTTTAGGCGCAAGCGGTCACTGGGATCCTAAAGGTGCAAAACAACATGGTTACCCAT"
    "GGCAAGATGATGCACACTTAGGTGATTTACCTGCATTAACTGTATTGCATGATGGCACAGCAACAAATCC"
    "TGTTTTAGCACCACGTCTTAAACATTTAGATGATGTTCGCGGTCACTCTATTATGATCCACACGGGTGGT"
    "GATAATCACTCCGATCATCCAGCTCCACTTGGCGGTGGCGGCCCACGTATGGCATGTGGCGTGATTAAAT"
    "A"
)
SODC_ACCESSION = "NC_013016.1:1293529-1294089"

_SOURCE_FILES = ("seroALLELES.fa", "ctrA.fa")

_TYPING = {
    "schema": "gapit.typing/2",
    "schemes": [
        {
            "name": "men_serogroup",
            "rules": [
                {
                    "model": "exact_set",
                    "phenotype": serogroup,
                    "requires": [f"SERO_{serogroup}"],
                    "excludes": [
                        f"SERO_{other}"
                        for other in ("A", "B", "C", "W", "X", "Y")
                        if other != serogroup
                    ],
                    "identity_floor": 90.0,
                    "coverage_floor": 90.0,
                    "notes": [],
                }
                for serogroup in ("A", "B", "C", "W", "X", "Y")
            ]
            + [
                {
                    "model": "exact_set",
                    "phenotype": "W/Y",
                    "requires": ["SERO_W", "SERO_Y"],
                    "excludes": ["SERO_A", "SERO_B", "SERO_C", "SERO_X"],
                    "identity_floor": 90.0,
                    "coverage_floor": 90.0,
                    "notes": [
                        "W/Y co-detection is not a final serogroup.",
                        "NEXT STEP [EX7E]: run menwy (MDU-PHL/meningotype) on the"
                        " assembly; it aligns synG and translates reference positions"
                        " 919-945 (synG-relative, not genome coordinates); P=W, G=Y,"
                        " S=W/Y. Do not report W/Y as final.",
                    ],
                }
            ],
            "cutoff": 1.0,
            "ambiguity_margin": 0.0,
            "fallback": "untypeable (atypical or non-meningitidis)",
            "control_gene": "sodC",
        }
    ],
}

_SCREEN = {
    "schema": "gapit.screen/1",
    "reads": {
        "min_breadth": 90.0,
        "min_identity": 80.0,
        "min_mapq": 20,
        "min_depth": 5.0,
        "min_relative_depth": 0.10,
        "relative_to": "sodC",
    },
    "contigs": {"minid": 80.0, "mincov": 80.0, "merge_fragments": True},
}

_FLOORS = {
    "schema": "gapit.floors/1",
    "default": 90.0,
    "genes": {"sodC": 80.0, "ctrA": 80.0},
}


def _records_for(fasta_path: Path, gene: str, function: str, product: str) -> Iterator[Record]:
    for record in iter_fasta(fasta_path):
        yield Record(
            db=NAME,
            gene=gene,
            sequence=record.sequence,
            accession=record.id,
            function=(function,),
            product=product.format(record.id),
            source_id=record.id,
        )


def transform(workdir: Path) -> Iterator[Record]:
    """Nucleotide records: serogroup amplicons, ctrA and sodC."""
    for record in iter_fasta(workdir / "seroALLELES.fa"):
        yield Record(
            db=NAME,
            gene=f"SERO_{record.id}",
            sequence=record.sequence,
            accession=record.id,
            function=("serotype",),
            product=f"Serogroup {record.id} capsule PCR amplicon",
            source_id=record.id,
        )
    yield from _records_for(workdir / "ctrA.fa", "ctrA", "control", "Capsule transport ctrA ({})")
    yield Record(
        db=NAME,
        gene="sodC",
        sequence=_SODC,
        accession=SODC_ACCESSION,
        function=("control",),
        product="Cu-Zn superoxide dismutase sodC (single-copy depth control)",
        source_id="sodC",
    )


def finalize(workdir: Path, build_dir: Path) -> None:
    """Write the serogroup typing, reads defaults and per-gene floors."""
    (build_dir / "typing.json").write_text(json.dumps(_TYPING, indent=2) + "\n", encoding="utf-8")
    (build_dir / "screen.json").write_text(json.dumps(_SCREEN, indent=2) + "\n", encoding="utf-8")
    (build_dir / "floors.json").write_text(json.dumps(_FLOORS, indent=2) + "\n", encoding="utf-8")


PROVIDER = Provider(
    name=NAME,
    description="Neisseria meningitidis serogroup + ctrA + sodC control (non-MLST)",
    vendor="MDU-PHL",
    source_urls=tuple(f"{RAW}/{name}" for name in _SOURCE_FILES),
    dbtype="nucl",
    transform=transform,
    finalize=finalize,
    license="GPL-3.0 (meningotype) + PubMLST Neisseria database terms",
    guidance=DbGuidance(
        references=(
            DbReference(
                citation=(
                    "Mothershed EA, Sacchi CT, Whitney AM, et al. Use of real-time PCR"
                    " to detect the capsule locus in Neisseria meningitidis."
                    " J Clin Microbiol 2004;42(1):320-328."
                ),
                doi="10.1128/JCM.42.1.320-328.2004",
            ),
            DbReference(
                citation=(
                    "Jolley KA, Brehony C, Maiden MC. Molecular typing of meningococci:"
                    " recommendations for target selection and nomenclature."
                    " FEMS Microbiol Rev 2007;31(1):89-96."
                ),
                doi="10.1111/j.1574-6976.2006.00057.x",
            ),
        ),
        caveats=(
            "Non-MLST only: use your gMLST pipeline for MLST/cc calls.",
            "Serogroup W/Y is called by co-detection of W and Y amplicons;"
            " EX7E follow-up is required for a final W vs Y call. Run menwy on the"
            " assembly; it aligns synG and translates reference positions 919-945"
            " (synG-relative, not genome coordinates): P=W, G=Y, S=W/Y.",
            "PorA/FetA/PorB/BAST allele-level finetyping is outside gapit's"
            " gene-presence scope; use the gMLST database for strain typing and"
            " an external finetyping pipeline for antigen allele IDs.",
            "sodC is a conserved single-copy Neisseria control; use its depth as a"
            " relative denominator, not as an absolute species-specific marker.",
            "Database content follows MDU-PHL/meningotype (GPL-3.0) and PubMLST"
            " Neisseria terms; fetch-on-demand only, not bundled.",
        ),
    ),
)
