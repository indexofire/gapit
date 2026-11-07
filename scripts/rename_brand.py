#!/usr/bin/env python3
"""One-shot brand rename gapit -> gspec for the gspec repository.

Protects every contract identifier that must not move in a 0.8.0 rename:

  * JSON-schema ids          gapit.<doc>/<n>
  * Python module paths      gapit.cli, gapit.errors, gapit.cli:app, ...
  * env vars                 GAPIT_DATADIR, GAPIT_ENGINE
  * on-disk names            gapit-manifest.json, the `gapit|` FASTA prefix,
                             the `gapit/v1` header-format id, the default
                             datadir `~/.local/share/gspec/db`
  * source/test identifiers  is_gapit_header, test_gapit_*, _gapit_

Only the standalone word `gapit` (prose, CLI name, version banner, vendor
strings) is rewritten.  Run with --check to see the diff size first.
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
    ".github",
}
SKIP_FILES = {
    "src/gapit/dbcodec.py",  # owns the `gapit|` prefix literal
    "scripts/rename_brand.py",
    "CHANGELOG.md",  # history keeps the name it shipped under
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
    ".cfg",
    ".ini",
    ".sh",
}

#: (pattern, placeholder) applied in order before the brand rewrite
PROTECT: list[tuple[str, str]] = [
    (r"\bgapit\.(?:[a-z_][a-z0-9_]*)/\d+", "SCHEMA"),  # gapit.typing/2
    (r"\bgapit\.(?:[a-z_][a-z0-9_]*)", "MODULE"),  # gapit.cli, gapit.errors
    (r"gapit-manifest(?:\.json)?", "MANIFEST"),
    (r"gapit\|", "HEADERPREFIX"),  # native FASTA prefix
    (r"gapit/v1", "HEADERFORMAT"),
    (r"share/gapit/", "DATADIRPATH"),
    (r"src/gapit/", "SRCDIR"),
    (r"gapit/data/", "DATADIR2"),
]


def rewrite(text: str) -> str:
    stash: list[str] = []

    def keep(match: re.Match[str]) -> str:
        stash.append(match.group(0))
        return f"\x00{len(stash) - 1}\x00"

    for pattern, _ in PROTECT:
        text = re.sub(pattern, keep, text)
    text = re.sub(r"(?<![A-Za-z0-9_])gapit(?![A-Za-z0-9_])", "gspec", text)
    for index, original in enumerate(stash):
        text = text.replace(f"\x00{index}\x00", original)
    return text


def targets() -> list[Path]:
    out: list[Path] = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix not in SUFFIXES:
            continue
        rel = path.relative_to(ROOT).as_posix()
        if any(part in SKIP_DIRS for part in path.relative_to(ROOT).parts):
            continue
        if rel in SKIP_FILES:
            continue
        out.append(path)
    return sorted(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="report only, write nothing")
    args = ap.parse_args()
    changed = []
    for path in targets():
        original = path.read_text(encoding="utf-8")
        updated = rewrite(original)
        if updated != original:
            changed.append(
                (path.relative_to(ROOT).as_posix(), original.count("gapit"), updated.count("gapit"))
            )
            if not args.check:
                path.write_text(updated, encoding="utf-8")
    for rel, before, after in changed:
        print(f"{rel:70s} gapit {before:4d} -> {after:4d}")
    print(f"\n{len(changed)} files {'would change' if args.check else 'changed'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
