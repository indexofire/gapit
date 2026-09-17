"""Summary parity harness: byte-diff `gapit summary` against real abricate
--summary over the committed synthetic report fixtures (SPEC.md §6).

Opt-in, like run_parity.py: run with `pixi run -e parity summary-parity`.
Uses only preexisting table fixtures under tests/data/summary — it never runs
screening and never fetches data.
"""

import difflib
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import NoReturn

PROJECT = Path(__file__).resolve().parents[2]
FIXTURES = Path(__file__).resolve().parent.parent / "data" / "summary"
DEFAULT_BIN = PROJECT / ".pixi" / "envs" / "default" / "bin"


def fail(message: str) -> NoReturn:
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def locate() -> tuple[Path, Path]:
    abricate = shutil.which("abricate")
    if abricate is None:
        fail("abricate not on PATH (run via: pixi run -e parity summary-parity)")
    gapit_bin = DEFAULT_BIN / "gapit"
    if not gapit_bin.is_file():
        fail(f"gapit CLI not found: {gapit_bin} (run: pixi install)")
    return Path(abricate), gapit_bin


def tool_env(abricate_bin: Path) -> dict[str, str]:
    env = os.environ.copy()
    env["PATH"] = os.pathsep.join([str(abricate_bin.parent), str(DEFAULT_BIN), env.get("PATH", "")])
    return env


def diff_case(
    label: str,
    abricate_argv: list[str],
    gapit_argv: list[str],
    env: dict[str, str],
    cwd: Path,
) -> bool:
    abricate_result = subprocess.run(
        ["abricate", *abricate_argv], cwd=cwd, env=env, capture_output=True, text=True, check=False
    )
    gapit_result = subprocess.run(
        [str(DEFAULT_BIN / "gapit"), *gapit_argv],
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    if abricate_result.returncode != 0 or gapit_result.returncode != 0:
        print(
            f"FAIL {label}: abricate rc={abricate_result.returncode}, "
            f"gapit rc={gapit_result.returncode}"
        )
        return False
    identical = abricate_result.stdout == gapit_result.stdout
    print(f"{'OK  ' if identical else 'DIFF'} {label}")
    if not identical:
        sys.stdout.writelines(
            difflib.unified_diff(
                abricate_result.stdout.splitlines(True),
                gapit_result.stdout.splitlines(True),
                fromfile=f"abricate/{label}",
                tofile=f"gapit/{label}",
            )
        )
    return identical


def main() -> int:
    abricate_bin, _gapit_bin = locate()
    env = tool_env(abricate_bin)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        cases: list[tuple[str, list[str], list[str], Path]] = [
            (
                "dutch",
                ["--summary", "multi_sample.tsv"],
                ["summary", "multi_sample.tsv"],
                FIXTURES,
            ),
            (
                "multi",
                ["--summary", "sample_a.tsv", "sample_b.tsv", "empty.tsv"],
                ["summary", "sample_a.tsv", "sample_b.tsv", "empty.tsv"],
                FIXTURES,
            ),
            (
                "identity",
                ["--summary", "--identity", "sample_a.tsv", "sample_b.tsv"],
                ["summary", "--identity", "--quiet", "sample_a.tsv", "sample_b.tsv"],
                FIXTURES,
            ),
            (
                "duplicate",
                ["--summary", "sample_a.tsv", "sample_a.tsv"],
                ["summary", "--quiet", "sample_a.tsv", "sample_a.tsv"],
                FIXTURES,
            ),
            (
                "csv",
                ["--summary", "--csv", "sample_a.csv", "sample_b.csv"],
                ["summary", "--format", "csv", "sample_a.csv", "sample_b.csv"],
                FIXTURES,
            ),
            (
                "nopath-key-sort",
                ["--summary", "--nopath", "1dir/zeta.tsv", "2dir/mid.tsv"],
                ["summary", "--nopath", "1dir/zeta.tsv", "2dir/mid.tsv"],
                work,
            ),
        ]
        # nopath-key-sort works on copies of the (empty) header-only fixture in
        # subdirectories, proving rows sort by full key while labels basename.
        header = (FIXTURES / "empty.tsv").read_text(encoding="utf-8")
        for name in ("1dir/zeta.tsv", "2dir/mid.tsv"):
            path = work / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(header, encoding="utf-8")
        failures = sum(
            0 if diff_case(label, abricate_argv, gapit_argv, env, cwd) else 1
            for label, abricate_argv, gapit_argv, cwd in cases
        )
    print(f"summary parity: {len(cases) - failures}/{len(cases)} cases byte-identical")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
