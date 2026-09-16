"""Regenerate the parity corpus from the real abricate databases.

Embeds exact, ~95%-identity mutated, and ~60% truncated copies of real genes
from abricate's ncbi/card databases, plus a junk contig. Deterministic: fixed
seed, fixed selection rule (first suitable genes in file order).

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


if __name__ == "__main__":
    main()
