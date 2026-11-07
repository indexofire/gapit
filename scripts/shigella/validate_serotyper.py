#!/usr/bin/env python3
"""Compare a gspec serotype typing table with an external reference.

Three reference shapes are auto-detected from the first line:

- **ShigaTyper** two-column TSV (``FILE<TAB>LABEL``, trailing whitespace
  tolerated);
- **ShigaPass** ``ShigaPass_summary.csv`` (semicolon-separated; composes the
  flexneri subtype from ``Predicted_FlexSerotype`` when the serotype is
  ``SF1-5``);
- **ShigeiFinder** results (``#SAMPLE``/``SEROTYPE`` TSV).

GapIt format: the normal ``gspec typing`` TSV (FILE, SCHEME, PHENOTYPE,
CONFIDENCE, ...).

Usage::

    python scripts/shigella/validate_serotyper.py \
        gapit_typing.tsv shigapass_summary.csv
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path

_ACCESSION = re.compile(r"^(GC[AF]_\d+)")
_SPECIES = ("boydii", "flexneri", "sonnei", "dysenteriae")
_LABEL = re.compile(
    r"^(?:S\.|Shigella)\s+(boydii|flexneri|sonnei|dysenteriae)(?:\s+(.*))?$",
    re.IGNORECASE,
)
_CODE = re.compile(r"^(SB|SD|SF|SS)\s*([0-9A-Za-z./()-]*)$", re.IGNORECASE)

#: Canonical tokens are lower-case with whitespace and parentheses removed.
_TOKEN_ALIASES = {
    "1c7a": "1c",
    "7a": "1c",
    "4c": "xv",
    "xv4c": "xv",
    "3batypicaloac1b": "3boac1b",
    "1rfp-": "1",
    "y": "y",
    "yv": "yv",
}
#: ShigaPass/ShigeiFinder numeric names for the ShigaTyper provisional types.
_NUMBER_ALIASES = {
    ("boydii", "22"): "provisionale1621-54",
    ("dysenteriae", "16"): "provisional96-265",
    ("dysenteriae", "18"): "provisionale670-74",
}
#: Reference labels that cannot be compared at the serotype level.
_BROAD = {"", "unknown", "none", "nd", "n/a", "1-5", "3/16"}


def _key(name: str) -> str:
    """A stable genome key across all three reference formats."""
    base = Path(name.strip()).name
    match = _ACCESSION.match(base)
    return match.group(1) if match is not None else base.removesuffix(".fna")


def _squash(text: str) -> str:
    return re.sub(r"[\s()]", "", text).lower()


def _matches(expected: tuple[str, str], got: tuple[str, str]) -> bool:
    """Exact serotype match, with a generic ``provisional`` reference
    accepting any specific provisional label for the same species."""
    if expected == got:
        return True
    return expected[1] == "provisional" and got[1].startswith("provisional")


def _canonical(label: str) -> tuple[str, str] | None:
    """(species, serotype token) or None when the label is undetermined."""
    text = " ".join(label.strip().split())
    low = text.lower()
    if not text or low in {"undetermined", "nd", "none", "n/a", "unknown", "eiec"}:
        return None
    if "not shigella" in low or low.startswith("shigella spp"):
        return None
    match = _LABEL.match(text)
    if match is not None:
        species = match.group(1).lower()
        rest = (match.group(2) or "").strip()
    else:
        code = _CODE.match(text)
        if code is None:
            return None
        prefix = code.group(1).upper()
        species = {"SB": "boydii", "SD": "dysenteriae", "SF": "flexneri", "SS": "sonnei"}[prefix]
        rest = code.group(2)
    token = _squash(rest)
    if "/" in token:
        # A reference that lists several serotypes is not comparable.
        return None
    token = _TOKEN_ALIASES.get(token, token)
    token = _NUMBER_ALIASES.get((species, token), token)
    if token in _BROAD:
        # Species-only labels are comparable for S. sonnei (the reference
        # does not distinguish form I/II) but not for the other species.
        return (species, "") if species == "sonnei" else None
    if species == "sonnei":
        return (species, "")
    return (species, token)


def _shigapass_label(serotype: str, flexneri: str) -> str:
    serotype = serotype.strip()
    flexneri = flexneri.strip()
    if serotype == "SF1-5":
        return f"S. flexneri {flexneri}" if flexneri and flexneri != "Unknown" else ""
    if serotype == "SF6":
        return "S. flexneri 6"
    if serotype == "SS":
        return "S. sonnei"
    if serotype.startswith("SB"):
        return f"S. boydii {serotype[2:]}"
    if serotype.startswith("SD"):
        return f"S. dysenteriae {serotype[2:]}"
    if serotype.startswith("SDprov"):
        return f"S. dysenteriae {serotype[2:]}"
    return serotype


def _load_shigapass(path: Path) -> dict[str, str]:
    reference: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter=";"):
            name = (row.get("Name") or "").strip()
            if not name:
                continue
            label = _shigapass_label(
                row.get("Predicted_Serotype") or "",
                row.get("Predicted_FlexSerotype") or "",
            )
            if label:
                reference[_key(name)] = label
    return reference


def _load_shigeifinder(path: Path) -> dict[str, str]:
    reference: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = (row.get("#SAMPLE") or "").strip()
            label = (row.get("SEROTYPE") or "").strip()
            if name and label:
                reference[_key(name)] = label
    return reference


def _load_two_column(path: Path) -> dict[str, str]:
    reference: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            continue
        reference[_key(fields[0])] = fields[1].strip()
    return reference


def _load_reference(path: Path) -> tuple[str, dict[str, str]]:
    first = path.read_text(encoding="utf-8").splitlines()[0]
    if first.startswith("Name;rfb;"):
        return "shigapass", _load_shigapass(path)
    if first.startswith("#SAMPLE\t"):
        return "shigeifinder", _load_shigeifinder(path)
    return "shigatyper", _load_two_column(path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gapit_typing", type=Path)
    parser.add_argument("reference", type=Path)
    args = parser.parse_args(argv)

    reference_format, reference = _load_reference(args.reference)
    ours: dict[str, tuple[str, str]] = {}
    with args.gapit_typing.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            ours[_key(row["FILE"])] = (row["PHENOTYPE"], row["CONFIDENCE"])

    pairs = [
        (name, reference[name], ours[name])
        for name in sorted(reference)
        if name in ours and _canonical(reference[name]) is not None
    ]
    counts: Counter[str] = Counter()
    mismatches: list[str] = []
    fallbacks: list[str] = []
    for name, expected, (got, confidence) in pairs:
        expected_key = _canonical(expected)
        assert expected_key is not None
        got_key = None if got == "-" else _canonical(got)
        if got == "-" or confidence == "ambiguous":
            counts["ambiguous"] += 1
        elif got_key is None:
            counts["fallback"] += 1
            fallbacks.append(f"{name}: expected {expected!r}, got {got!r}")
        elif _matches(expected_key, got_key):
            counts["exact"] += 1
        elif got_key[1] == "6/10":
            counts["ambiguous"] += 1
        else:
            counts["mismatch"] += 1
            mismatches.append(f"{name}: expected {expected!r}, got {got!r}")
    total = len(pairs)
    print(f"reference: {reference_format}")
    print(f"determined reference genomes: {total}")
    for key in ("exact", "ambiguous", "fallback", "mismatch"):
        value = counts[key]
        print(f"  {key:10s} {value:4d}  {100 * value / total if total else 0:.1f}%")
    if fallbacks:
        print("fallbacks:")
        print("\n".join(f"  {line}" for line in fallbacks))
    if mismatches:
        print("mismatches:")
        print("\n".join(f"  {line}" for line in mismatches))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
