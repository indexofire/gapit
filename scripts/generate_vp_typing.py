#!/usr/bin/env python3
"""Generate the committed vp_ok typing asset (gapit.typing/2) from the GBKs.

A DEVELOPER tool, not a CLI surface: recomputes the discriminating-gene
analysis over the vpautils OAgc + CPSgc GenBank files (nothing is hardcoded
from the development matrix — the matrix only documents the method) and
emits the three-scheme document installed by the vp_ok cluster bundle
(o_group + k_group + serotype compose). Support helpers (homology pass,
name pairing, JSON/split writers) live in scripts/vp_typing_support.py.

Method (per typed locus OL*/KL*; OLU*/KLU* get no rule — the OUT/KUT
fallbacks own them):

- genes are keyed by their build-time RESOLVED ids (resolve_gene_ids,
  exactly what a gapit.cluster/2 build runs), so key-gene references
  validate against the built features;
- tier-1 key genes = the locus's sequence-unique genes (>=95% identity
  homology nowhere else; support module), antigen-processing genes first
  (wzx/wzy/wzm/wzt/wcaA/wza/wzb/wzc/kpsM/kpsT), then longest CDS, capped
  at 4;
- HARD loci (no sequence-unique gene; the development matrix flags
  OL1/OL1V1/OL1V2/OL4V1/OL10/OL11/KL71, trusted over the recomputation via
  --matrix) key on the locus's antigen genes and rely on the 90%
  coverage / 95% identity floors against homologous siblings; THIN loci (one
  unique gene) key on that gene. Every caveat lands in the rule notes.

Variants fold into their base antigen (OL4V1 -> O4, KL20V1 -> K20); the
vpautils reference defines 12 O groups (no OL13), so O3 stays plain O3.

Usage:
    python scripts/generate_vp_typing.py OAgc.gbk CPSgc.gbk \
        --output src/gapit/data/typing/vp_ok.json \
        --matrix scripts/vp_rules_matrix.tsv \
        [--calibration-dir DIR]

The committed scripts/vp_rules_matrix.tsv carries the development-matrix
HARD flags (loci whose every CDS has a close cross-locus homolog) — the
only matrix signal the generator consumes.
"""

import argparse
import re
import sys
import tempfile
from pathlib import Path

from vp_typing_support import (
    Gene,
    gene_feature_names,
    matrix_hard_loci,
    non_unique_genes,
    write_json,
    write_splits,
)

from gapit.clusterbuild import parse_cluster_input, resolve_gene_ids
from gapit.gbfeatures import LocusSequence
from gapit.typing_models import TypingDocument, TypingScheme
from gapit.typing_rules import (
    ClusterMatchRule,
    CoverageComponent,
    IdentityComponent,
    KeyGenesComponent,
)

ANTIGEN_GENES = ("wzx", "wzy", "wzm", "wzt", "wcaA", "wza", "wzb", "wzc", "kpsM", "kpsT")
FAMILY_KEYS = {
    "O": ("wzx", "wzy", "wzm", "wzt", "wcaA"),
    "K": ("wza", "wzb", "wzc", "kpsM", "kpsT"),
}
LOCUS_ID = re.compile(r"^([OK])L(\d+)(?:V(\d+))?$")
MAX_KEY_GENES = 4
HARD_NOTE = (
    "HARD: no sequence-unique CDS (sibling loci are close homologs); key genes are the"
    " locus's antigen-processing genes, present when the locus is truly carried —"
    " homologous siblings rely on the 90% coverage / 95% identity floors"
)
THIN_NOTE = "THIN: single sequence-unique CDS; key genes are that gene alone"

# Calibration-derived rule overrides (real-genome harness, batch_test corpus):
# every entry fixes a DEMONSTRATED misfire and lands in the rule notes.
#   OL5 — divergent O5 lineages (GCA_000732985/736325, locus 100% covered at
#     95.3% identity) carry waaA/degT alleles below the gene identity floor;
#     keys retuned to O-locus backbone genes present in those carriers.
#   KL12, KL15 — partial-deletion carriers sit at ~88% locus breadth
#     (GCA_015808295.2 KL12 88.4%, GCA_000707465.2 KL15 88.9%); coverage
#     floor relaxed to 85 and keys retuned to genes present in them.
# Deliberately NOT overridden (documented limitations, see the calibration
# report): KL28 (carrier at 77.5% breadth — declined, vpautils calls it
# Medium), OL12/KLU59 (GCA_041681985.1 carries a 51.7%-breadth truncated
# O locus and an untyped KLU; OUT/KUT is the honest answer).
CALIBRATION_OVERRIDES: dict[str, dict[str, object]] = {
    "OL5": {
        "keys": ["OL5_004", "OL5_008", "OL5_011", "OL5_027"],
        "note": (
            "calibration: keys retuned to backbone genes (rfaF/rfbA/rffA/rfaD) — divergent"
            " O5 lineages carry waaA/degT alleles below the 90% gene identity floor"
            " (GCA_000732985.1/GCA_000736325.1)"
        ),
    },
    "KL12": {
        "coverage_floor": 85.0,
        "keys": ["KL12_002", "KL12_031", "KL12_027", "KL12_025"],
        "note": (
            "calibration: coverage floor 85 (partial-deletion carriers at ~88% breadth,"
            " GCA_015808295.2) and keys retuned to genes present in them"
        ),
    },
    "KL15": {
        "coverage_floor": 85.0,
        "keys": ["KL15_010", "KL15_011", "KL15_006", "KL15_028"],
        "note": (
            "calibration: coverage floor 85 (partial-deletion carrier at ~89% breadth,"
            " GCA_000707465.2) and keys retuned to genes present in it"
        ),
    },
}


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _sort_key(locus: str) -> tuple[str, int, int, str]:
    match = LOCUS_ID.fullmatch(locus)
    assert match is not None, locus
    return (match.group(1), int(match.group(2)), int(match.group(3) or 0), locus)


def load_genes(*gbks: Path) -> tuple[list[Gene], dict[str, str]]:
    """Every typed locus's genes in O-then-K file order, with resolved ids
    (the cluster/2 build transform) and the raw-id map for split docs."""
    names = gene_feature_names(*gbks)
    loci: list[LocusSequence] = []
    for gbk in gbks:
        loci.extend(parse_cluster_input(gbk))
    genes: list[Gene] = []
    raw_of_resolved: dict[str, str] = {}
    for locus, fixed in zip(loci, resolve_gene_ids(loci), strict=True):
        match = LOCUS_ID.fullmatch(locus.features.id)
        if match is None:
            continue
        for raw, gene in zip(locus.features.genes, fixed.features.genes, strict=True):
            piece = locus.sequence[gene.start - 1 : gene.end]
            genes.append(
                Gene(
                    locus=locus.features.id,
                    family=match.group(1),
                    raw_id=raw.gene_id,
                    resolved_id=gene.gene_id,
                    name=names.get((locus.features.id, gene.locus_tag or ""), raw.gene_id),
                    sequence=piece if gene.strand == "+" else _reverse_complement(piece),
                )
            )
            raw_of_resolved[gene.gene_id] = raw.gene_id
    return genes, raw_of_resolved


def key_genes(locus_genes: list[Gene], shared: set[str]) -> tuple[list[str], str]:
    """One typed locus's tier-1 key genes (resolved ids) plus caveat note:
    antigen genes first then longest CDS over the unique set; HARD loci fall
    back to the locus's own antigen genes; THIN loci key on their one gene."""
    unique = [gene for gene in locus_genes if gene.resolved_id not in shared]
    if len(unique) == 1:
        return [unique[0].resolved_id], THIN_NOTE
    if unique:
        ordered = sorted(
            unique,
            key=lambda gene: (
                gene.name not in ANTIGEN_GENES,
                -len(gene.sequence),
                gene.resolved_id,
            ),
        )
        return [gene.resolved_id for gene in ordered[:MAX_KEY_GENES]], ""
    preference = FAMILY_KEYS[locus_genes[0].family]
    antigen_carrying = sorted(
        (gene for gene in locus_genes if gene.name in ANTIGEN_GENES),
        key=lambda gene: (
            preference.index(gene.name) if gene.name in preference else len(preference),
            gene.resolved_id,
        ),
    )
    if antigen_carrying:
        return [gene.resolved_id for gene in antigen_carrying[:MAX_KEY_GENES]], HARD_NOTE
    longest = sorted(locus_genes, key=lambda gene: (-len(gene.sequence), gene.resolved_id))
    return [gene.resolved_id for gene in longest[:MAX_KEY_GENES]], (
        HARD_NOTE + " (no named antigen gene: longest CDS)"
    )


def cluster_rule(
    locus: str, keys: list[str], note: str, coverage_floor: float = 90.0
) -> ClusterMatchRule:
    phenotype = LOCUS_ID.fullmatch(locus)
    notes: list[str] = []
    if phenotype.group(3):
        notes.append(f"variant locus {locus} folds into {phenotype.group(1)}{phenotype.group(2)}")
    if f"{phenotype.group(1)}{phenotype.group(2)}" == "O3":
        notes.append("vpautils reference defines 12 O groups (no O13 locus); plain O3")
    if note:
        notes.append(note)
    return ClusterMatchRule(
        model="cluster_match",
        phenotype=f"{phenotype.group(1)}{phenotype.group(2)}",
        coverage=CoverageComponent(weight=0.55, floor=coverage_floor),
        identity=IdentityComponent(weight=0.25, floor=95.0),
        key_genes=KeyGenesComponent(weight=0.2, genes=keys),
        notes=notes,
    )


def build_document(
    genes: list[Gene], shared: set[str], matrix_hard: set[str]
) -> tuple[TypingDocument, dict[str, list[ClusterMatchRule]]]:
    """The three-scheme gapit.typing/2 document, rules sorted by family,
    antigen number, then variant ordinal; also the per-scheme rule lists for
    the calibration splits."""
    by_locus: dict[str, list[Gene]] = {}
    for gene in genes:
        by_locus.setdefault(gene.locus, []).append(gene)
    rules: dict[str, list[ClusterMatchRule]] = {"O": [], "K": []}
    for locus in sorted(by_locus, key=_sort_key):
        match = LOCUS_ID.fullmatch(locus)
        recomputed_hard = all(gene.resolved_id in shared for gene in by_locus[locus])
        locus_shared = (
            shared
            if locus not in matrix_hard and not recomputed_hard
            else {gene.resolved_id for gene in by_locus[locus]}
        )
        override = CALIBRATION_OVERRIDES.get(locus, {})
        override_keys = override.get("keys")
        if isinstance(override_keys, list):
            keys, note = list(override_keys), str(override.get("note", ""))
        else:
            keys, note = key_genes(by_locus[locus], locus_shared)
        rules[match.group(1)].append(
            cluster_rule(locus, keys, note, float(override.get("coverage_floor", 90.0)))
        )
    schemes = [
        TypingScheme(
            name=f"{family.lower()}_group",
            rules=tuple(rules[family]),
            cutoff=0.9,
            ambiguity_margin=0.05,
            fallback="OUT" if family == "O" else "KUT",
        )
        for family in ("O", "K")
    ]
    schemes.append(
        TypingScheme(
            name="serotype",
            compose="{o_group}:{k_group}",
            cutoff=0.9,
            ambiguity_margin=0.05,
            fallback="OUT:KUT",
        )
    )
    document = TypingDocument.model_validate({"schemes": [s.model_dump() for s in schemes]})
    return document, rules


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="generate the vp_ok typing asset")
    parser.add_argument("o_gbk", type=Path, help="OAgc.gbk (vpautils serotype database)")
    parser.add_argument("k_gbk", type=Path, help="CPSgc.gbk (vpautils serotype database)")
    parser.add_argument("--output", type=Path, required=True, help="merged vp_ok.json path")
    parser.add_argument("--matrix", type=Path, default=None, help="development matrix (HARD flags)")
    parser.add_argument("--calibration-dir", type=Path, default=None, help="split-doc output dir")
    args = parser.parse_args(argv)

    genes, raw_of_resolved = load_genes(args.o_gbk, args.k_gbk)
    with tempfile.TemporaryDirectory(prefix="vp-typing.") as workdir:
        shared = non_unique_genes(genes, Path(workdir))
    document, rules = build_document(genes, shared, matrix_hard_loci(args.matrix))
    write_json(document, args.output)
    if args.calibration_dir is not None:
        write_splits(rules, raw_of_resolved, args.calibration_dir)
    hard = sorted(
        rule.phenotype
        for family_rules in rules.values()
        for rule in family_rules
        if any("HARD" in note for note in rule.notes)
    )
    typed_loci = {gene.locus for gene in genes}
    print(
        f"wrote {args.output}: {sum(len(v) for v in rules.values())} rules over"
        f" {len(typed_loci)} typed loci (HARD: {','.join(hard) or 'none'})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
