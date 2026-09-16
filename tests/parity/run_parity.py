"""Parity harness: byte-diff gapit against real abricate on the committed corpus.

Run with: pixi run -e parity parity
Requires the parity environment (abricate) and the default environment (gapit);
the harness runs on the default env's python and merges both envs' PATHs so
both tools share the same blastn binary.
"""

import difflib
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import NoReturn

from gapit.db import make_blast_db

PROJECT = Path(__file__).resolve().parents[2]
CORPUS = Path(__file__).resolve().parent / "corpus"
DEFAULT_BIN = PROJECT / ".pixi" / "envs" / "default" / "bin"
DB_NAMES = ("ncbi", "card")


def fail(message: str) -> NoReturn:
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def locate() -> tuple[Path, Path, Path]:
    """(abricate binary, abricate datadir, gapit binary)."""
    abricate = shutil.which("abricate")
    if abricate is None:
        fail("abricate not on PATH (run via: pixi run -e parity parity)")
    abricate_bin = Path(abricate).resolve()
    datadir = abricate_bin.parents[1] / "db"
    if not datadir.is_dir():
        fail(f"abricate datadir not found: {datadir}")
    gapit_bin = DEFAULT_BIN / "gapit"
    if not gapit_bin.is_file():
        fail(f"gapit CLI not found: {gapit_bin} (run: pixi install)")
    return abricate_bin, datadir, gapit_bin


def tool_env(abricate_bin: Path) -> dict[str, str]:
    """PATH with the parity env first (abricate's perl must win over the perl
    that ships with blast in the default env); both tools then share the parity
    env's blastn/any2fasta, keeping the comparison same-binary."""
    env = os.environ.copy()
    env["PATH"] = os.pathsep.join([str(abricate_bin.parent), str(DEFAULT_BIN), env.get("PATH", "")])
    return env


def ensure_indices(datadir: Path) -> None:
    for name in DB_NAMES:
        db_dir = datadir / name
        has_index = (db_dir / "sequences.nin").exists() or (db_dir / "sequences.pin").exists()
        if not has_index:
            print(f"building missing BLAST index for {name} ...")
            make_blast_db(db_dir / "sequences", name)


def run_tool(argv: list[str], env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(argv, cwd=CORPUS, env=env, capture_output=True, text=True, check=False)


def main() -> int:
    abricate_bin, datadir, gapit_bin = locate()
    ensure_indices(datadir)
    env = tool_env(abricate_bin)
    corpora = sorted(CORPUS.glob("*.fa"))
    if not corpora:
        fail(f"no corpus files in {CORPUS} (run tests/parity/make_corpus.py)")
    failures = 0
    for db_name in DB_NAMES:
        for corpus in corpora:
            abricate_result = run_tool(
                ["abricate", "--db", db_name, "--datadir", str(datadir), "--nopath", corpus.name],
                env,
            )
            gapit_result = run_tool(
                [
                    str(gapit_bin),
                    "screen",
                    "--db",
                    db_name,
                    "--datadir",
                    str(datadir),
                    "--nopath",
                    corpus.name,
                ],
                env,
            )
            if abricate_result.returncode != 0 or gapit_result.returncode != 0:
                failures += 1
                print(
                    f"FAIL {db_name} {corpus.name}: "
                    f"abricate rc={abricate_result.returncode}, gapit rc={gapit_result.returncode}"
                )
                continue
            identical = abricate_result.stdout == gapit_result.stdout
            print(f"{'OK  ' if identical else 'DIFF'} {db_name} {corpus.name}")
            if not identical:
                failures += 1
                sys.stdout.writelines(
                    difflib.unified_diff(
                        abricate_result.stdout.splitlines(True),
                        gapit_result.stdout.splitlines(True),
                        fromfile=f"abricate/{db_name}/{corpus.name}",
                        tofile=f"gapit/{db_name}/{corpus.name}",
                    )
                )
    total = len(DB_NAMES) * len(corpora)
    print(f"parity: {total - failures}/{total} cases byte-identical")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
