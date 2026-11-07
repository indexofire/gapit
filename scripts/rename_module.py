#!/usr/bin/env python3
"""One-shot module rename gapit -> gspec (0.8.1) for the gspec repository.

Why this exists: 0.8.0 shipped the code under ``gapit/``, so a metadata-only
``gapit`` transition package could not work -- pip installs ``gspec`` first,
then uninstalls the old ``gapit``, and that uninstall deletes every path the
old RECORD lists, including gspec's freshly installed ``gapit/*`` files.

0.8.1 moves the import package to ``gspec/`` and leaves the ``gapit`` import to
the transition package itself.

Protects the frozen contract spellings: ``gapit.<doc>/<n>`` schema ids,
``gapit-manifest.json``, the ``gapit|`` FASTA prefix, ``gapit/v1``, the legacy
``GAPIT_*`` env vars and identifiers such as ``is_gapit_header``.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path.home() / "repos/my-project/gapit"
SKIP_DIRS = {
    ".git",
    ".pixi",
    "site",
    "__pycache__",
    ".ruff_cache",
    ".pytest_cache",
    "dist",
    "packaging",
}
SKIP_FILES = {
    "scripts/rename_module.py",
    "scripts/rename_brand.py",
    "CHANGELOG.md",
}
SUFFIXES = {
    ".py",
    ".md",
    ".toml",
    ".yml",
    ".yaml",
    ".json",
    ".jsonl",
    ".cff",
    ".txt",
    ".sh",
    ".cfg",
    ".ini",
}

PROTECT = [
    r"\bgapit\.(?:[a-z_][a-z0-9_]*)/\d+",  # gapit.typing/2 (schema id)
    r"gapit-manifest(?:\.json)?",
    r"gapit\|",
    r"gapit/v1",
    r"GAPIT_[A-Z_]+",
]

REWRITES = [
    (r"(?<![\w.])gapit\.", "gspec."),  # module attribute / dotted path
    (r"\bfrom gapit import\b", "from gspec import"),
    (r"\bimport gapit\b", "import gspec"),
    (r"(?<![\w/-])gapit(?=[/\"' ,)])", "gspec"),
]


def rewrite(text: str) -> str:
    stash: list[str] = []

    def keep(match: re.Match[str]) -> str:
        stash.append(match.group(0))
        return f"\x00{len(stash) - 1}\x00"

    for pattern in PROTECT:
        text = re.sub(pattern, keep, text)
    for pattern, repl in REWRITES:
        text = re.sub(pattern, repl, text)
    for index, original in enumerate(stash):
        text = text.replace(f"\x00{index}\x00", original)
    return text


def targets() -> list[Path]:
    out = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix not in SUFFIXES:
            continue
        rel = path.relative_to(ROOT)
        if any(part in SKIP_DIRS for part in rel.parts):
            continue
        if rel.as_posix() in SKIP_FILES:
            continue
        out.append(path)
    return sorted(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    changed = 0
    for path in targets():
        original = path.read_text(encoding="utf-8")
        updated = rewrite(original)
        if updated == original:
            continue
        changed += 1
        if not args.check:
            path.write_text(updated, encoding="utf-8")
    print(f"{changed} files {'would change' if args.check else 'changed'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
