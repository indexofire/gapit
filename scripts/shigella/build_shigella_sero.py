#!/usr/bin/env python3
"""Build the bundled Shigella marker panel from ShigaTyper's public-domain FASTA.

This is the data-curation step, not a serotype resolver.  It rewrites the
ShigaTyper marker names into gspec/v1 headers, emits marker metadata for the
future profile engine, and writes the bundle manifest/screen defaults.

Usage::

    python scripts/shigella/build_shigella_sero.py \
        /path/to/ShigellaRef5.fasta \
        --dest src/gspec/data/dbs/shi_sero

The source file is ShigaTyper 2.0.5's ``ShigellaRef5.fasta`` (Public Domain).
ShigaPass GPL content is deliberately not copied; only the ShigaTyper
public-domain marker panel is used in this first slice.
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import UTC, datetime
from pathlib import Path

from gspec.dbcodec import encode_seqid
from gspec.fasta import iter_fasta

_NAME = re.compile(r"^(?P<family>S[bd f s]|S[bd f s]|Sd|Sb|Sf|Ss)(?P<num>\d+)?(?P<rest>.*)$", re.I)


def _category(marker: str) -> str:
    lower = marker.lower()
    if lower in {"ipah_c", "ipab"} or lower in {"eclacy", "cada"}:
        return "virulence" if lower.startswith("ipa") else "biochemical"
    if any(token in lower for token in ("wzx", "wzy", "rfp", "wbam", "heparinase")):
        return "o_antigen"
    if any(token in lower for token in ("gtr", "oac", "xv", "opt")):
        return "poac"
    if "methylase" in lower:
        return "o_antigen_modification"
    return "other"


def _serotype_group(marker: str) -> str | None:
    lower = marker.lower()
    if lower in {"ipah_c", "ipab", "eclacy", "cada"}:
        return None
    match = re.match(r"^(sb|sd|sf|ss)(\d+)?", lower)
    if match is None:
        return None
    family, number = match.group(1).upper(), match.group(2)
    if family == "SF":
        return "SF1-5" if marker in {"Sf_wzx", "Sf_wzy"} else "SF6"
    if family == "SS":
        return "SS"
    if number is None:
        return family + "Prov"
    return f"{family}{number}"


def _product(marker: str) -> str:
    return f"Shigella serotype marker {marker}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="ShigaTyper ShigellaRef5.fasta")
    parser.add_argument("--dest", type=Path, required=True)
    parser.add_argument("--name", default="shi_sero")
    args = parser.parse_args(argv)

    records = list(iter_fasta(args.source))
    if not records:
        raise SystemExit(f"no FASTA records in {args.source}")
    if not any(record.id == "ipaH_c" for record in records):
        raise SystemExit(f"{args.source} has no ipaH_c marker")
    dest: Path = args.dest
    dest.mkdir(parents=True, exist_ok=True)
    markers: list[dict[str, object]] = []
    with (dest / "sequences").open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            marker = record.id
            # ShigaTyper's ipaH_c is one allele; divergent Shigella ipaH
            # alleles retain different parts of it.  Split the control into
            # two non-overlapping presence fragments: keeping the full-length
            # record would lose a fragment to the screen's per-query dedup
            # (the full-length hit claims the same contig interval, then
            # fails the coverage filter).  Both fragments keep the gene name
            # `ipaH_c`, so the existing control gate folds them.
            parts: list[tuple[str, str, str | None]] = (
                [
                    ("ipaH_c", record.sequence[:640], "ipaH_c:1-640"),
                    ("ipaH_c", record.sequence[640:780], "ipaH_c:641-780"),
                ]
                if marker == "ipaH_c"
                else [(marker, record.sequence, None)]
            )
            for gene, sequence, derived_from in parts:
                accession = (
                    f"ShigaTyper:ShigellaRef5:{derived_from[7:]}"
                    if derived_from is not None
                    else "ShigaTyper:ShigellaRef5"
                )
                header = encode_seqid(args.name, gene, accession, ("serotype",))
                handle.write(f">{header}\n")
                for offset in range(0, len(sequence), 60):
                    handle.write(sequence[offset : offset + 60] + "\n")
                entry: dict[str, object] = {
                    "gene": gene,
                    "category": _category(gene),
                    "serotype_group": _serotype_group(gene),
                    "length": len(sequence),
                }
                if derived_from is not None:
                    entry["derived_from"] = derived_from
                markers.append(entry)
    (dest / "markers.json").write_text(
        json.dumps(
            {
                "schema": "gspec.shigella_markers/1",
                "source": {
                    "tool": "ShigaTyper",
                    "version": "2.0.5",
                    "file": "ShigellaRef5.fasta",
                    "license": "Public Domain",
                },
                "markers": markers,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (dest / "bundled.json").write_text(
        json.dumps(
            {
                "name": args.name,
                "description": (
                    "Shigella molecular serotype marker panel (ShigaTyper public-domain markers)"
                ),
                "vendor": "gspec-curated (ShigaTyper public-domain markers)",
                "dbtype": "nucl",
                "kind": "gene",
                "snapshotted": datetime.now(UTC).strftime("%Y-%m-%d"),
                "license": "Public Domain (ShigaTyper marker sequences)",
                "note": (
                    "ShigaTyper 2.0.5 ShigellaRef5 marker panel, reformatted as gspec/v1;"
                    " markers.json carries the O-antigen/POAC/virulence roles and"
                    " serotype.json carries the first conservative marker-set resolver."
                    " The ipaH_c control is split into two non-overlapping presence"
                    " fragments (1-640 and 641-780), a gspec-derived sensitivity"
                    " marker for divergent ipaH alleles."
                ),
                "guidance": {
                    "references": [
                        {
                            "citation": (
                                "Wu Y, Lau HK, Lee T, Lau DK, Payne J. In Silico"
                                " Serotyping Based on Whole-Genome Sequencing Improves"
                                " the Accuracy of Shigella Identification. Appl Environ"
                                " Microbiol. 2019;85(7):e00165-19."
                            ),
                            "url": "https://github.com/CFSAN-Biostatistics/shigatyper",
                            "note": "ShigaTyper marker panel; Public Domain.",
                        },
                        {
                            "citation": (
                                "Senchenkova SN, Feng L, Yang J, et al. Structural and"
                                " Genetic Characterization of the Shigella boydii Type 10"
                                " and Type 6 O Antigens. J Bacteriol."
                                " 2005;187(7):2551-2554."
                            ),
                            "doi": "10.1128/JB.187.7.2551-2554.2005",
                            "note": (
                                "wbaM is intact in type 10 and interrupted by IS629"
                                " in type 6; the 5' alignment start separates them."
                            ),
                        },
                    ],
                    "caveats": [
                        (
                            "Contig screens resolve S. boydii 6 vs 10 from the"
                            " wbaM 5' alignment start; reads tables carry no"
                            " alignment geometry and leave that pair unresolved."
                        ),
                        (
                            "The ipaH_c control is two non-overlapping"
                            " presence fragments (1-640 and 641-780), not the"
                            " full-length ShigaTyper gene: divergent ipaH"
                            " alleles retain different parts of the gene."
                        ),
                        (
                            "ShigaPass GPL3 database content is deliberately not"
                            " copied into this MIT bundle."
                        ),
                    ],
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (dest / "screen.json").write_text(
        json.dumps(
            {
                "schema": "gspec.screen/1",
                "contigs": {"minid": 80.0, "mincov": 50.0, "merge_fragments": True},
                "reads": {"min_breadth": 50.0, "min_identity": 80.0},
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"wrote {dest}: {len(markers)} markers")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
