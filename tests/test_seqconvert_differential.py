"""Differential tests: native seqconvert vs the real any2fasta binary.

Compares PARSED records (id, description, sequence) — not raw bytes: gapit
re-wraps sequence at 60 columns while the perl keeps input wrapping, a
difference blast cannot see. Skipped unless any2fasta is on PATH; the difftest
pixi env provides it:

    pixi run -e difftest difftest
"""

import shutil
import subprocess
from pathlib import Path

import pytest

from gapit.seqconvert import to_fasta_lines

CONVERT = Path(__file__).parent / "data" / "convert"
ANY2FASTA = shutil.which("any2fasta")
FIXTURES = [
    "sample.fa",
    "sample.fq",
    "sample.gbk",
    "sample.embl",
    "sample.fa.gz",
    "sample.gbk.gz",
    "sample.embl.bz2",
]

pytestmark = pytest.mark.skipif(
    ANY2FASTA is None, reason="any2fasta not on PATH (run via: pixi run -e difftest difftest)"
)


def parsed(text: str) -> list[tuple[str, str, str]]:
    """(id, description, sequence) triples from FASTA text."""
    records: list[tuple[str, str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append(_split(header, "".join(chunks)))
            header, chunks = line[1:], []
        else:
            chunks.append(line)
    if header is not None:
        records.append(_split(header, "".join(chunks)))
    return records


def _split(header: str, sequence: str) -> tuple[str, str, str]:
    parts = header.split(None, 1)
    return (parts[0] if parts else "", parts[1] if len(parts) > 1 else "", sequence)


def run_any2fasta(path: Path) -> str:
    assert ANY2FASTA is not None  # narrowed for the checker; pytestmark skips otherwise
    result = subprocess.run(
        [ANY2FASTA, "-q", "-u", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


@pytest.mark.parametrize("name", FIXTURES)
def test_matches_any2fasta_parsed_records(name: str) -> None:
    """Given a convert fixture, When converted natively and by the real
    `any2fasta -q -u`, Then both yield identical parsed records (id,
    description, uppercase sequence) for every record."""
    fixture = CONVERT / name
    expected = parsed(run_any2fasta(fixture))
    native = parsed("".join(to_fasta_lines(fixture)))
    assert native == expected
    assert expected, "fixture produced no records"
