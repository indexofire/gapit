#!/usr/bin/env python3
"""Generate the first gapit.serotype/1 rule set for the ShigaTyper marker panel.

This is a conservative marker-set port: it covers the unambiguous wzx/wzy
markers, the S. sonnei forms, the S. flexneri gtr/Oac/Xv combinations, the
Sd1 rfp variant, and the S. boydii 6 vs 10 wbaM mechanism. Cases that need
evidence the assembly table does not carry stay explicitly ambiguous.

Usage::

    python scripts/shigella/build_serotype_rules.py \
        src/gapit/data/dbs/shi_sero
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from gapit.serotype_models import SerotypeDocument, SerotypeRule

POAC = ("gtrI", "gtrII", "gtrIV", "gtrV", "gtrX", "gtrIC", "Oac", "Oac1b", "Xv")
#: IS629 interrupts S. boydii serotype 6 wbaM at ~252/957 (26.3%); 20% leaves
#: margin for assembly soft-clipping while still requiring the type-10 5' end.
WBA_SB6_START_PCT = 20.0
FLEX = {
    "S. flexneri 1a": ("gtrI",),
    "S. flexneri 1b": ("gtrI", "Oac1b"),
    "S. flexneri 1c (7a)": ("gtrI", "gtrIC"),
    "S. flexneri 7b": ("gtrI", "gtrIC", "Oac1b"),
    "S. flexneri 2a": ("gtrII",),
    "S. flexneri 2av": ("gtrII", "Xv"),
    "S. flexneri 2b": ("gtrII", "gtrX"),
    "S. flexneri 3a": ("gtrX", "Oac"),
    "S. flexneri 3b": ("Oac",),
    "S. flexneri 3b (oac1b)": ("Oac1b",),
    "S. flexneri 4a": ("gtrIV",),
    "S. flexneri 4av": ("gtrIV", "Xv"),
    "S. flexneri 4b": ("gtrIV", "Oac"),
    "S. flexneri 4bv": ("gtrIV", "Oac", "Xv"),
    "S. flexneri 5a": ("gtrV", "Oac"),
    "S. flexneri 5b": ("gtrV", "gtrX", "Oac"),
    "S. flexneri X": ("gtrX",),
    "S. flexneri Xv (4c)": ("gtrX", "Xv"),
    "S. flexneri Yv": ("Xv",),
}


def _rule(phenotype: str, **kwargs: object) -> SerotypeRule:
    return SerotypeRule.model_validate({"phenotype": phenotype, **kwargs})


def _rules() -> list[SerotypeRule]:
    rules: list[SerotypeRule] = []
    # S. dysenteriae 1: rfp status is a real subcall.
    rules.append(
        _rule(
            "S. dysenteriae 1",
            requires=("Sd1_wzx", "Sd1_rfp"),
            notes=("ShigaTyper: Sd1 rfp+",),
        )
    )
    rules.append(
        _rule(
            "S. dysenteriae 1 (rfp-)",
            requires=("Sd1_wzx",),
            excludes=("Sd1_rfp",),
            notes=("ShigaTyper: Sd1 wzx+ / rfp- (phenotypically negative)",),
        )
    )
    # S. boydii 1 vs 20 (galF/heparinase marker).
    rules.append(
        _rule(
            "S. boydii 20",
            requires=("Sb1_wzx", "heparinase"),
            notes=("ShigaTyper: Sb1 rfb + heparinase = SB20",),
        )
    )
    rules.append(
        _rule(
            "S. boydii 1",
            requires=("Sb1_wzx",),
            excludes=("heparinase",),
            notes=("ShigaTyper: Sb1 rfb without heparinase",),
        )
    )
    # S. sonnei forms.
    rules.append(
        _rule(
            "S. sonnei form I",
            requires=("Ss_wzx", "Ss_methylase"),
            notes=(
                "ShigaTyper: Ss_wzx + Ss_methylase + cadA; cadA is optional"
                " because S. sonnei isolates can lose the cadA call while"
                " retaining the form I O-antigen genes",
            ),
        )
    )
    rules.append(
        _rule(
            "S. sonnei form II",
            requires=("Ss_methylase", "cadA"),
            excludes=("Ss_wzx",),
            notes=("ShigaTyper: cadA + Ss_methylase without Ss_wzx",),
        )
    )
    # S. flexneri 6 is a separate O backbone.
    rules.append(
        _rule(
            "S. flexneri 6",
            requires=("Sf6_wzx",),
            notes=("ShigaTyper: Sf6_wzx",),
        )
    )
    # S. flexneri 1-5/X/Y subserotypes: exact POAC set on the shared Sf rfb.
    for phenotype, required in FLEX.items():
        rules.append(
            _rule(
                phenotype,
                requires=("Sf_wzx", *required),
                excludes=tuple(gene for gene in POAC if gene not in required),
                notes=("ShigaTyper POAC set",),
            )
        )
    rules.append(
        _rule(
            "S. flexneri Y",
            requires=("Sf_wzx",),
            excludes=POAC,
            notes=("ShigaTyper: Sf rfb with no POAC modification",),
        )
    )
    # Direct Sd/Sb wzx markers. Multiple matching rules become ambiguous.
    for number in range(2, 16):
        rules.append(_rule(f"S. dysenteriae {number}", requires=(f"Sd{number}_wzx",)))
    rules.append(_rule("S. dysenteriae provisional 96-265", requires=("SdProv_wzx",)))
    rules.append(_rule("S. dysenteriae provisional E670-74", requires=("SdProvE_wzx",)))
    for number in range(2, 20):
        if number in (6, 10, 13):
            continue
        rules.append(_rule(f"S. boydii {number}", requires=(f"Sb{number}_wzx",)))
    rules.append(
        _rule(
            "S. boydii 10",
            requires=("Sb6_wzx", "wbaM"),
            max_alignment_start_pct={"wbaM": WBA_SB6_START_PCT},
            notes=(
                "Senchenkova et al. 2005: type 10 has an intact wbaM "
                "(ribofuranosyltransferase); a wbaM alignment covering the 5' end "
                "separates it from type 6",
            ),
        )
    )
    rules.append(
        _rule(
            "S. boydii 6",
            requires=("Sb6_wzx", "wbaM"),
            max_coverage={"wbaM": 95.0},
            min_alignment_start_pct={"wbaM": WBA_SB6_START_PCT},
            notes=(
                "Senchenkova et al. 2005: type 6 wbaM is interrupted by IS629; "
                "a 5'-truncated wbaM alignment (IS629 insertion at ~252/957) identifies type 6",
            ),
        )
    )
    rules.append(
        _rule(
            "S. boydii 6/10",
            requires=("Sb6_wzx",),
            excludes=("wbaM",),
            notes=("wbaM absent/weak; S. boydii 6 vs 10 cannot be separated",),
        )
    )
    rules.append(
        _rule(
            "S. boydii 13",
            requires=("Sb13_wzx",),
            control_exempt=True,
            notes=("SB13 is no longer considered Shigella; control-exempt by ShigaTyper",),
        )
    )
    rules.append(_rule("S. boydii provisional E1621-54", requires=("SbProv_wzx",)))
    return rules


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("db_dir", type=Path)
    args = parser.parse_args(argv)
    db_dir: Path = args.db_dir
    markers = json.loads((db_dir / "markers.json").read_text(encoding="utf-8"))["markers"]
    genes = {entry["gene"] for entry in markers}
    document = SerotypeDocument(
        name="shi_sero",
        control_gene="ipaH_c",
        fallback="Not Shigella/EIEC",
        rules=tuple(_rules()),
    )
    for rule in document.rules:
        for gene in (*rule.requires, *rule.requires_any, *rule.excludes):
            if gene not in genes:
                raise SystemExit(f"rule {rule.phenotype!r} references missing marker {gene!r}")
    (db_dir / "serotype.json").write_text(
        document.model_dump_json(indent=2, by_alias=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {db_dir / 'serotype.json'}: {len(document.rules)} rules")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
