#!/usr/bin/env python3
"""Support helpers of the vp_ok typing generator (analysis side).

Everything the generator (scripts/generate_vp_typing.py) needs that is not
document assembly: the cross-locus homology pass, the ``gene``-feature name
pairing, the deterministic JSON writer, the calibration split documents,
and the development matrix's HARD-flag reader.

- the homology pass marks a gene SHARED when some gene of ANOTHER locus
  aligns to it at >=95% identity over >=95% of its length (one all-vs-all
  minimap2 asm20 pass, single-threaded for reproducibility);
- the Bakta vpautils files name antigen genes only on ``gene`` features,
  never on the CDS (whose ids the cluster parser derives from
  /locus_tag), so names are paired onto CDS by locus_tag;
- the calibration harness (scripts/cluster_calibration.py) refuses
  multi-scheme databases, so per-scheme tuning runs over SPLIT databases:
  one single-scheme gapit.typing/2 document per family, re-keyed on
  PRIMARY gene ids (an unsplit build never resolves ids).
"""

import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

from gapit.paf import parse_paf_row
from gapit.typing_models import TypingDocument
from gapit.typing_rules import ClusterMatchRule, KeyGenesComponent

GENE_FEATURE = re.compile(r"^ {5}gene\s+\S.*?(?=^ {5}\S|\Z)", re.M | re.S)
HOMOLOG_IDENTITY = 0.95
HOMOLOG_COVERAGE = 0.95


@dataclass(frozen=True, slots=True)
class Gene:
    """One typed locus CDS: locus, family, primary and resolved ids, the
    ``gene``-feature name (antigen vocabulary), and the strand-oriented
    nucleotide sequence."""

    locus: str
    family: str
    raw_id: str
    resolved_id: str
    name: str
    sequence: str


def gene_feature_names(*gbks: Path) -> dict[tuple[str, str], str]:
    """(locus id, locus_tag) -> /gene name from ``gene`` features."""
    names: dict[tuple[str, str], str] = {}
    for gbk in gbks:
        for record in gbk.read_text(encoding="utf-8").split("LOCUS       ")[1:]:
            locus = record.split(maxsplit=1)[0]
            for block in GENE_FEATURE.findall(record):
                tags = re.findall(r'/locus_tag="([^"]+)"', block)
                genes = re.findall(r'/gene="([^"]+)"', block)
                if tags and genes:
                    names[(locus, tags[0])] = genes[0]
    return names


def non_unique_genes(genes: list[Gene], workdir: Path) -> set[str]:
    """Resolved ids with a >=95%-identity, >=95%-coverage cross-locus
    homolog (the deterministic all-vs-all minimap2 asm20 pass)."""
    fasta = workdir / "genes.fna"
    with fasta.open("w", encoding="utf-8", newline="\n") as handle:
        for gene in genes:
            handle.write(f">{gene.resolved_id}\n{gene.sequence}\n")
    result = subprocess.run(
        ["minimap2", "-x", "asm20", "-t", "1", "-N", "20", "-p", "0", str(fasta), str(fasta)],
        capture_output=True,
        text=True,
        check=True,
    )
    locus_of = {gene.resolved_id: gene.locus for gene in genes}
    shared: set[str] = set()
    for line in result.stdout.splitlines():
        if not line:
            continue
        row = parse_paf_row(line)
        if row.qname == row.tname or locus_of[row.qname] == locus_of[row.tname]:
            continue
        identity = row.nmatch / row.alen if row.alen else 0.0
        coverage = (row.qend - row.qstart) / row.qlen if row.qlen else 0.0
        if identity >= HOMOLOG_IDENTITY and coverage >= HOMOLOG_COVERAGE:
            shared.add(row.qname)
    return shared


def write_json(document: TypingDocument, path: Path) -> None:
    """Deterministic asset bytes: 2-space JSON + trailing newline, with
    empty notes lists and compose rules arrays stripped for readability."""
    data = json.loads(document.model_dump_json(by_alias=True))
    for scheme in data["schemes"]:
        if scheme.get("compose") is not None:
            del scheme["rules"]
        for rule in scheme.get("rules", ()):
            if not rule["notes"]:
                del rule["notes"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_splits(
    rules: dict[str, list[ClusterMatchRule]], raw_of_resolved: dict[str, str], outdir: Path
) -> None:
    """vp_ok_o.json / vp_ok_k.json beside the merged asset: the per-family
    rule lists re-keyed on primary gene ids for the calibration harness."""
    outdir.mkdir(parents=True, exist_ok=True)
    for family, fallback in (("O", "OUT"), ("K", "KUT")):
        swapped = [
            rule.model_copy(
                update={
                    "key_genes": KeyGenesComponent(
                        weight=rule.key_genes.weight,
                        genes=[raw_of_resolved[gene] for gene in rule.key_genes.genes],
                    )
                }
            )
            for rule in rules[family]
        ]
        split = TypingDocument.model_validate(
            {
                "schemes": [
                    {
                        "name": f"{family.lower()}_group",
                        "rules": [rule.model_dump() for rule in swapped],
                        "cutoff": 0.9,
                        "ambiguity_margin": 0.05,
                        "fallback": fallback,
                    }
                ]
            }
        )
        write_json(split, outdir / f"vp_ok_{family.lower()}.json")


def matrix_hard_loci(matrix: Path | None) -> set[str]:
    """Locus ids the development matrix flags HARD (notes column)."""
    if matrix is None:
        return set()
    return {
        line.split("\t", 2)[1]
        for line in matrix.read_text(encoding="utf-8").splitlines()[1:]
        if re.search(r"\bHARD", line)
    }
