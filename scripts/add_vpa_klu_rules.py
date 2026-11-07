#!/usr/bin/env python3
"""Append locus-only KUT rules for the updated vpa_ok KLU loci.

The main generator (``scripts/generate_vp_typing.py``) only emits rules for
typed ``OL*``/``KL*`` loci.  The current vpa_ok snapshot also carries
``KLU61``-``KLU100``: Kaptive treats these as typeable unknown loci, while
gspec represents the same state as ``K_GROUP=KUT`` plus the detected locus
and coverage.  This deterministic post-step installs one ``cluster_match``
rule per KLU locus (phenotype ``KUT``, all locus genes as key genes) so the
fallback evidence and ``unknown_locus`` status are available.

Usage::

    python scripts/add_vpa_klu_rules.py \
        [src/gapit/data/dbs/vpa_ok/typing.json] \
        [src/gapit/data/dbs/vpa_ok/features.json]

The operation is idempotent: a KLU locus already covered by any K-scheme rule
is left untouched.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from gapit.gbfeatures import FeaturesDocument
from gapit.typing_models import read_typing_document, validate_references


def _parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "typing",
        nargs="?",
        type=Path,
        default=root / "src/gapit/data/dbs/vpa_ok/typing.json",
    )
    parser.add_argument(
        "features",
        nargs="?",
        type=Path,
        default=root / "src/gapit/data/dbs/vpa_ok/features.json",
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    typing_path: Path = args.typing
    features_path: Path = args.features
    document = json.loads(typing_path.read_text(encoding="utf-8"))
    features = FeaturesDocument.model_validate_json(features_path.read_text(encoding="utf-8"))
    k_scheme = next(scheme for scheme in document["schemes"] if scheme["name"] == "k_group")
    covered = {
        gene.split("_", 1)[0] for rule in k_scheme["rules"] for gene in rule["key_genes"]["genes"]
    }
    klu_loci = sorted(
        (
            (int(locus.label[3:]), locus)
            for locus in features.loci
            if locus.label.startswith("KLU") and locus.label[3:].isdigit()
        ),
        key=lambda item: item[0],
    )
    added = 0
    for number, locus in klu_loci:
        label = locus.label
        if not 61 <= number <= 100 or label in covered:
            continue
        k_scheme["rules"].append(
            {
                "model": "cluster_match",
                "phenotype": "KUT",
                "coverage": {"weight": 0.55, "floor": 90.0},
                "identity": {"weight": 0.25, "floor": 95.0},
                "key_genes": {
                    "weight": 0.2,
                    "genes": [gene.gene_id for gene in locus.genes],
                },
                "notes": [f"locus-only: {label} has no known K phenotype; KUT with locus evidence"],
            }
        )
        added += 1
    if added:
        typing_path.write_text(
            json.dumps(document, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
    parsed = read_typing_document(typing_path)
    validate_references(parsed, features)
    print(f"{typing_path}: added {added} KLU rules; validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
