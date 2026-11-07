#!/usr/bin/env python3
"""Build the bundled V. parahaemolyticus VPaI-1..VPaI-10 cluster database.

``vpa_pai`` is a gspec-curated cluster bundle: ten pathogenicity-island
loci, each one island, with the CDS features used for gene-level evidence.
Sources are the public INSDC records behind Boyd et al. 2008
(PMID 18590559, BMC Microbiol 8:110):

- VPaI-1..VPaI-7 from RIMD 2210633 (NC_004603.1 / NC_004605.1);
- VPaI-8 and VPaI-9 from AQ3810 (GCA_000154045.1 scaffolds DS180230 and
  DS180088);
- VPaI-10 is the T6SS island (VP1387-VP1420) in RIMD 2210633.  The original
  paper labels the VPaI-9/VPaI-10 insertion sites in its Fig. S4 but does not
  give a separate VPaI-10 coordinate table; this boundary follows the T6SS
  region and is recorded as a caveat in bundled.json.

Usage::

    python scripts/vibrio/build_vpa_pai.py \
        --rimd RIMD2210633.gbff --aq3810 AQ3810.gbff \
        --dest src/gspec/data/dbs/vpa_pai
"""

from __future__ import annotations

import argparse
import json
import re
from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path

from gspec.gbfeatures import GeneFeature, LocusFeatures, write_features_document

#: (locus id, source record, start, end).  Coordinates are 1-based and
#: inclusive.  VPaI-1..7 use the published RIMD coordinates (the 2026
#: re-analysis table agrees with the 2008 ORF ranges); VPaI-8/9 use the
#: AQ3810 scaffold ORF spans; VPaI-10 is the T6SS CDS span.
LOCI: tuple[tuple[str, str, int, int], ...] = (
    ("VPaI-1", "NC_004603", 381054, 403433),
    ("VPaI-2", "NC_004603", 664402, 674289),
    ("VPaI-3", "NC_004603", 1121252, 1152668),
    ("VPaI-4", "NC_004603", 2240007, 2256166),
    ("VPaI-5", "NC_004603", 3084846, 3099979),
    ("VPaI-6", "NC_004605", 1325821, 1352643),
    ("VPaI-7", "NC_004605", 1390967, 1501509),
    ("VPaI-8", "DS180230", 6191, 22985),
    ("VPaI-9", "DS180088", 4122, 22794),
    ("VPaI-10", "NC_004603", 1471101, 1527506),
)

_LOCATION = re.compile(r"^(?P<complement>complement\()?(?P<rest>.*?)(?P<close>\))?$")


def _strip_quotes(value: str) -> str:
    value = value.strip()
    return value[1:-1] if len(value) >= 2 and value[0] == value[-1] == '"' else value


def parse_gbk(path: Path) -> dict[str, tuple[str, list[dict[str, object]]]]:
    """``record id -> (sequence, CDS records)`` for the simple locations we
    need.  A CDS record keeps its raw location and qualifier map."""
    records: dict[str, tuple[str, list[dict[str, object]]]] = {}
    record_id: str | None = None
    sequence: list[str] = []
    cds: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    section = ""
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.startswith("LOCUS"):
            if record_id is not None:
                records[record_id] = ("".join(sequence).upper(), cds)
            record_id = raw.split()[1]
            sequence = []
            cds = []
            current = None
            section = ""
            continue
        if record_id is None:
            continue
        if raw.startswith("ORIGIN"):
            section = "origin"
            continue
        if raw.startswith("//"):
            records[record_id] = ("".join(sequence).upper(), cds)
            record_id = None
            continue
        if section == "origin":
            sequence.append(re.sub(r"[^ACGTNacgtn]", "", raw[10:]))
            continue
        if raw.startswith("FEATURES"):
            section = "features"
            continue
        if section != "features":
            continue
        if raw.startswith("     CDS"):
            current = {"location": raw[21:].strip(), "quals": {}}
            cds.append(current)
        elif current is not None:
            stripped = raw.strip()
            if stripped.startswith("/"):
                key, _, value = stripped[1:].partition("=")
                current["quals"].setdefault(key, []).append(_strip_quotes(value))  # type: ignore[union-attr]
            elif stripped and not current["quals"]:  # type: ignore[union-attr]
                current["location"] += stripped  # type: ignore[union-attr]
            elif stripped:
                key = next(reversed(current["quals"]))  # type: ignore[union-attr]
                current["quals"][key][-1] += stripped  # type: ignore[union-attr]
    return records


def location_span(location: str) -> tuple[int, int, str]:
    """Min/max 1-based span and strand for simple/complement/join locations."""
    strand = "-" if location.startswith("complement(") else "+"
    inner = location
    if inner.startswith("complement(") and inner.endswith(")"):
        inner = inner[11:-1]
    parts = inner[5:-1].split(",") if inner.startswith("join(") and inner.endswith(")") else [inner]
    spans = []
    for part in parts:
        match = re.fullmatch(r"[<>]?(\d+)\.\.[<>]?(\d+)", part.strip())
        if match is None:
            raise ValueError(f"unsupported CDS location: {location!r}")
        spans.append((int(match.group(1)), int(match.group(2))))
    return min(a for a, _ in spans), max(b for _, b in spans), strand


def _qual(quals: dict[str, list[str]], key: str) -> str:
    return (quals.get(key) or [""])[0].strip()


def _gene_id(quals: dict[str, list[str]], fallback: str) -> str:
    for key in ("gene", "old_locus_tag", "locus_tag"):
        value = _qual(quals, key)
        if value:
            return value
    return fallback


def build_locus(
    locus_id: str,
    sequence: str,
    cds: list[dict[str, object]],
    start: int,
    end: int,
) -> tuple[LocusFeatures, str]:
    """One locus slice plus its fully-contained CDS, relative coordinates."""
    locus_sequence = sequence[start - 1 : end]
    genes: list[GeneFeature] = []
    used: set[str] = set()
    for index, feature in enumerate(cds, start=1):
        location = str(feature["location"])
        try:
            gene_start, gene_end, strand = location_span(location)
        except ValueError:
            continue
        if gene_start < start or gene_end > end:
            continue
        quals = feature["quals"]  # type: ignore[assignment]
        gene_id = _gene_id(quals, f"{locus_id}_{index:03d}")
        if gene_id in used:
            suffix = 2
            while f"{gene_id}_{suffix}" in used:
                suffix += 1
            gene_id = f"{gene_id}_{suffix}"
        used.add(gene_id)
        genes.append(
            GeneFeature(
                gene_id=gene_id,
                start=gene_start - start + 1,
                end=gene_end - start + 1,
                strand="-" if strand == "-" else "+",
                product=_qual(quals, "product"),
            )
        )
    return (
        LocusFeatures(id=locus_id, label=locus_id, type="PAI", genes=tuple(genes)),
        locus_sequence,
    )


def _write_fasta(loci: Iterable[tuple[str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for locus_id, sequence in loci:
            handle.write(f">{locus_id}\n")
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset : offset + 60] + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rimd", type=Path, required=True)
    parser.add_argument("--aq3810", type=Path, required=True)
    parser.add_argument("--dest", type=Path, required=True)
    args = parser.parse_args(argv)

    genomes: dict[str, tuple[str, list[dict[str, object]]]] = {}
    genomes.update(parse_gbk(args.rimd))
    genomes.update(parse_gbk(args.aq3810))

    features: list[LocusFeatures] = []
    sequences: list[tuple[str, str]] = []
    for locus_id, record_id, start, end in LOCI:
        if record_id not in genomes:
            raise SystemExit(f"record {record_id!r} for {locus_id} not found")
        sequence, cds = genomes[record_id]
        if end > len(sequence):
            raise SystemExit(f"{locus_id} end exceeds {record_id} length")
        locus_features, locus_sequence = build_locus(locus_id, sequence, cds, start, end)
        features.append(locus_features)
        sequences.append((locus_id, locus_sequence))

    dest: Path = args.dest
    dest.mkdir(parents=True, exist_ok=True)
    write_features_document(tuple(features), dest / "features.json")
    _write_fasta(sequences, dest / "sequences")
    (dest / "bundled.json").write_text(
        json.dumps(
            {
                "name": "vpa_pai",
                "description": "V. parahaemolyticus VPaI-1..VPaI-10 pathogenicity islands",
                "vendor": "gspec-curated (Boyd et al. 2008; INSDC)",
                "dbtype": "nucl",
                "kind": "cluster",
                "snapshotted": datetime.now(UTC).strftime("%Y-%m-%d"),
                "license": "Public domain (INSDC sequences)",
                "note": (
                    "Ten VPaI loci reconstructed from RIMD 2210633"
                    " (NC_004603.1/NC_004605.1) and AQ3810"
                    " (GCA_000154045.1); VPaI-10 follows the RIMD T6SS region"
                    " VP1387-VP1420 (VP1386 is absent from the RefSeq annotation)."
                    " Built as one locus per island through the"
                    " cluster pipeline; no typing scheme."
                ),
                "guidance": {
                    "references": [
                        {
                            "citation": (
                                "Boyd EF, Cohen ALV, Naughton LM, Ussery DW,"
                                " Binnewies TT, Stine OC, Parent MA. Molecular"
                                " analysis of the emergence of pandemic Vibrio"
                                " parahaemolyticus. BMC Microbiol. 2008;8:110."
                            ),
                            "doi": "10.1186/1471-2180-8-110",
                            "note": (
                                "Defines VPaI-1..VPaI-10, including the AQ3810"
                                " islands VPaI-8/VPaI-9 and the VPaI-9/VPaI-10"
                                " insertion-site comparison (Fig. S4)."
                            ),
                        },
                        {
                            "citation": (
                                "Hurley CC, Quirke A, Reen FJ, Boyd EF. Four"
                                " genomic islands that mark post-1995 pandemic"
                                " Vibrio parahaemolyticus isolates. BMC Genomics."
                                " 2006;7:104."
                            ),
                            "doi": "10.1186/1471-2164-7-104",
                            "note": "The original VPaI-1..VPaI-7 island description.",
                        },
                    ],
                    "caveats": (
                        "VPaI-1..7 coordinates follow RIMD 2210633; VPaI-8/9"
                        " follow the AQ3810 ORF descriptions. VPaI-10 is the"
                        " RIMD T6SS region VP1387-VP1420 (VP1386 is absent from the"
                        " RefSeq annotation): the 2008 paper shows"
                        " the VPaI-9/VPaI-10 insertion-site comparison but does"
                        " not publish a separate VPaI-10 coordinate table, so"
                        " this boundary is the documented T6SS reconstruction.",
                        "Island presence and coverage are not a pathogenicity"
                        " verdict: VPaIs can be partial, rearranged, or lost, and"
                        " several islands occur in pre-1995 and environmental"
                        " strains. Read the per-gene evidence, not just the locus"
                        " coverage.",
                        "The panel carries no typing scheme; it answers which"
                        " VPaIs are present and how completely, not which"
                        " pathotype or serotype the isolate is.",
                    ),
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"wrote {dest}: {len(features)} loci")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
