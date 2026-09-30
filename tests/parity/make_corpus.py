"""Regenerate the parity corpus from the real abricate databases.

Embeds exact, ~95%-identity mutated, and ~60% truncated copies of real genes
from abricate's ncbi/card databases, plus a junk contig. Deterministic: fixed
seed, fixed selection rule (first suitable genes in file order).

Corpus files (ncbi_a/ncbi_b = first two suitable ncbi genes in file order,
card_a = first suitable card gene; suitable = unique name, 300..900 nt):

- 01_exact_and_junk.fa: exact ncbi_a + card_a copies, plus a junk contig.
- 02_mutated.fa: every-20th-base mutants (~95% identity) of ncbi_a and card_a.
- 03_truncated.fa: exact ncbi_b plus its first 60% (coverage 60% -> filtered).
- 04_reverse.fa: reverse complements (case-preserving) of the same genes —
  exact revcomp of ncbi_a and card_a (minus-strand hits: STRAND '-', swapped
  subject coords) and an every-20th-base mutant of ncbi_a, revcomped
  (minus strand + mismatch positions in the COVERAGE_MAP).
- 05_threshold.fa: boundaries around the default --minid 80/--mincov 80 —
  ncbi_b truncated to 79% of its length (coverage just under mincov ->
  filtered by both tools) and to 81% (kept), plus a seeded-random 21%
  divergence copy of card_a (~79% identity, under blastn -perc_identity ->
  no hit row at all). Divergence placement is random rather than periodic so
  exact runs long enough to seed blastn (-task blastn, word 11) survive and
  the identity filter, not seed absence, is what drops the hit.

Run with the default env's python (requires the parity env to be installed):
    .pixi/envs/default/bin/python tests/parity/make_corpus.py
"""

import random
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
DB_DIR = PROJECT / ".pixi" / "envs" / "parity" / "db"
CORPUS = Path(__file__).resolve().parent / "corpus"
MIN_LEN = 300
MAX_LEN = 900

FLIP = {"A": "C", "C": "G", "G": "T", "T": "A", "a": "c", "c": "g", "g": "t", "t": "a"}


def read_genes(path: Path) -> list[tuple[str, str]]:
    """(gene, sequence) pairs from an abricate-format sequences file."""
    genes: list[tuple[str, str]] = []
    gene = ""
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if gene:
                genes.append((gene, "".join(chunks)))
            header = line[1:]
            gene = header.split("~~~")[1] if "~~~" in header else header
            chunks = []
        else:
            chunks.append(line.strip())
    if gene:
        genes.append((gene, "".join(chunks)))
    return genes


def pick_genes(genes: list[tuple[str, str]], count: int) -> list[tuple[str, str]]:
    """First `count` genes within the size window, in file order, skipping
    repeated gene names (databases contain duplicate-name entries)."""
    picked: list[tuple[str, str]] = []
    seen: set[str] = set()
    for gene, sequence in genes:
        if MIN_LEN <= len(sequence) <= MAX_LEN and gene not in seen:
            picked.append((gene, sequence))
            seen.add(gene)
        if len(picked) == count:
            break
    return picked


def mutate(sequence: str, every: int) -> str:
    """Substitute every `every`-th base (5% divergence at every=20)."""
    letters = list(sequence)
    for i in range(every // 2, len(letters), every):
        letters[i] = FLIP.get(letters[i], letters[i])
    return "".join(letters)


def revcomp(sequence: str) -> str:
    """Reverse complement, preserving case (minus-strand embedding)."""
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "a": "t", "c": "g", "g": "c", "t": "a"}
    return "".join(complement.get(base, base) for base in reversed(sequence))


def diverge(sequence: str, rate: float, seed: int) -> str:
    """Seed-random substitution of a `rate` fraction of the bases. Unlike the
    periodic mutate(), placement clumps, leaving exact runs long enough to
    seed blastn — so the -perc_identity filter, not seed absence, drops the
    resulting alignment."""
    rng = random.Random(seed)
    letters = list(sequence)
    for i in sorted(rng.sample(range(len(letters)), int(rate * len(letters)))):
        letters[i] = FLIP.get(letters[i], letters[i])
    return "".join(letters)


def junk(length: int) -> str:
    return "".join(random.Random(20260916).choices("ACGT", k=length))


def contig_id(prefix: str, gene: str) -> str:
    return prefix + "".join(char if char.isalnum() else "_" for char in gene)


def write_fasta(path: Path, contigs: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for seqid, sequence in contigs:
            handle.write(f">{seqid}\n")
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i : i + 60] + "\n")
    print(f"wrote {path} ({[len(seq) for _, seq in contigs]} nt)")


def main() -> None:
    ncbi = pick_genes(read_genes(DB_DIR / "ncbi" / "sequences"), 2)
    card = pick_genes(read_genes(DB_DIR / "card" / "sequences"), 1)
    if len(ncbi) < 2 or len(card) < 1:
        raise SystemExit("not enough suitable genes found in abricate databases")
    ncbi_a, ncbi_b = ncbi
    (card_a,) = card
    CORPUS.mkdir(parents=True, exist_ok=True)
    write_fasta(
        CORPUS / "01_exact_and_junk.fa",
        [
            (contig_id("exact_ncbi_", ncbi_a[0]), ncbi_a[1]),
            (contig_id("exact_card_", card_a[0]), card_a[1]),
            ("junk_contig", junk(250)),
        ],
    )
    write_fasta(
        CORPUS / "02_mutated.fa",
        [
            (contig_id("mut95_ncbi_", ncbi_a[0]), mutate(ncbi_a[1], 20)),
            (contig_id("mut95_card_", card_a[0]), mutate(card_a[1], 20)),
        ],
    )
    write_fasta(
        CORPUS / "03_truncated.fa",
        [
            (contig_id("exact_ncbi_", ncbi_b[0]), ncbi_b[1]),
            (contig_id("trunc60_ncbi_", ncbi_b[0]), ncbi_b[1][: int(0.6 * len(ncbi_b[1]))]),
        ],
    )
    write_fasta(
        CORPUS / "04_reverse.fa",
        [
            (contig_id("revcomp_ncbi_", ncbi_a[0]), revcomp(ncbi_a[1])),
            (contig_id("revcomp_card_", card_a[0]), revcomp(card_a[1])),
            (contig_id("revcomp95_ncbi_", ncbi_a[0]), revcomp(mutate(ncbi_a[1], 20))),
        ],
    )
    write_fasta(
        CORPUS / "05_threshold.fa",
        [
            (contig_id("cov79_ncbi_", ncbi_b[0]), ncbi_b[1][: int(0.79 * len(ncbi_b[1]))]),
            (contig_id("cov81_ncbi_", ncbi_b[0]), ncbi_b[1][: int(0.81 * len(ncbi_b[1]))]),
            (contig_id("mut79_card_", card_a[0]), diverge(card_a[1], 0.21, 20260917)),
        ],
    )


if __name__ == "__main__":
    main()
