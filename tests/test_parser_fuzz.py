"""Deterministic parser fuzz smoke: malformed boundaries stay typed.

No hypothesis dependency: a fixed-seed byte/syntax mutator exercises the
BLAST/PAF/cs/GBK/GFF parsers with 100 variants each. Every outcome must be a
typed GapitError (or a successful parse), never a raw IndexError/KeyError/
ValueError/ValidationError from a user-controlled input.
"""

import random
from collections.abc import Callable
from pathlib import Path

from gapit.blast import parse_blast_row
from gapit.cluster_math import parse_cs
from gapit.errors import GapitError
from gapit.gbfeatures import parse_genbank_features
from gapit.gffparse import parse_gff3_features
from gapit.paf import parse_paf_row

GBK = """LOCUS       TEST1                    60 bp    DNA     linear   BCT 01-JAN-2020
FEATURES             Location/Qualifiers
     CDS             1..30
                     /gene="abc"
ORIGIN
        1 acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt
//
"""

GFF = """##gff-version 3
L1\tx\tCDS\t1\t30\t.\t+\t0\tgene=abc
##FASTA
>L1
acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt
"""

BLAST = "\t".join(
    [
        "ctg",
        "1",
        "10",
        "100",
        "id",
        "1",
        "10",
        "100",
        "plus",
        "1e-20",
        "10",
        "100.0",
        "0",
        "0",
        "prod",
    ]
)
PAF = "\t".join(
    [
        "read",
        "100",
        "0",
        "100",
        "+",
        "target",
        "100",
        "0",
        "100",
        "100",
        "100",
        "60",
        "tp:A:P",
        "NM:i:0",
        "cs:Z::100",
    ]
)
CS = ":10*AC-ACG+TT:5"

CHARS = ' \t\n\r<>()..,;/="abc0123456789LOCUSORIGINFEATURESCDS##FASTA:+-*'


def _mutate(rng: random.Random, text: str) -> str:
    chars = list(text)
    for _ in range(rng.randint(1, 10)):
        op = rng.choice(("del", "ins", "rep"))
        position = rng.randrange(len(chars)) if chars else 0
        if op == "del" and chars:
            chars.pop(position)
        elif op == "ins":
            chars.insert(position, rng.choice(CHARS))
        elif chars:
            chars[position] = rng.choice(CHARS)
    return "".join(chars)


def _text_cases() -> tuple[tuple[str, str, Callable[[str], object]], ...]:
    return (
        ("blast", BLAST, parse_blast_row),
        ("paf", PAF, parse_paf_row),
        ("cs", CS, lambda value: parse_cs(value, 0)),
    )


def _file_cases() -> tuple[tuple[str, str, str, Callable[[Path], object]], ...]:
    return (
        ("gbk", GBK, ".gbk", parse_genbank_features),
        ("gff", GFF, ".gff3", parse_gff3_features),
    )


def test_mutated_text_parsers_never_raise_raw_exceptions() -> None:
    for name, text, parser in _text_cases():
        rng = random.Random(f"gapit-parser-fuzz-{name}")
        for _ in range(100):
            try:
                parser(_mutate(rng, text))
            except GapitError:
                pass
            except Exception as exc:  # pragma: no cover - failure detail
                raise AssertionError(f"{name} raised raw {type(exc).__name__}: {exc}") from exc


def test_mutated_file_parsers_never_raise_raw_exceptions(tmp_path: Path) -> None:
    for name, text, suffix, parser in _file_cases():
        rng = random.Random(f"gapit-parser-fuzz-{name}")
        for _ in range(100):
            path = tmp_path / f"input{suffix}"
            path.write_text(_mutate(rng, text), encoding="utf-8")
            try:
                parser(path)
            except GapitError:
                pass
            except Exception as exc:  # pragma: no cover - failure detail
                raise AssertionError(f"{name} raised raw {type(exc).__name__}: {exc}") from exc
