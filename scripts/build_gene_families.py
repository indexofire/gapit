#!/usr/bin/env python3
"""Build the optional panel gene-family sidecar for a cluster database.

The expensive part is one all-vs-all protein search over the panel's genes;
it is a database-build step, never a per-sample cost.  The resulting
``gene_families.json`` drives containment diagnostics and stable family
colours in the locus map.

Usage::

    python scripts/build_gene_families.py DB_DIR [--output PATH]
        [--threads N] [--min-identity 90] [--min-coverage 80]

The script expects ``features.json`` and ``sequences`` in ``DB_DIR`` and
requires ``makeblastdb`` + ``blastp`` on PATH (BLAST+ is already a gspec
dependency).
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import time
from collections.abc import Iterator
from pathlib import Path

from gspec.fasta import iter_fasta
from gspec.gbfeatures import FeaturesDocument
from gspec.gene_families import (
    GENE_FAMILY_FILENAME,
    HomologyHit,
    build_family_document,
    cluster_families,
    filter_homology_hits,
    write_gene_families,
)
from gspec.proctools import run_tool

# Standard genetic code (bacterial table 11 uses the same internal codons;
# the alternative start codons only affect the first residue).
_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODON = {
    first + second + third: _AMINO[index]
    for index, (first, second, third) in enumerate(
        (a, b, c) for a in _BASES for b in _BASES for c in _BASES
    )
}
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def _translate(sequence: str) -> str:
    """Translate a CDS from the first base, stopping at the first stop."""
    amino = []
    for index in range(0, len(sequence) - 2, 3):
        residue = _CODON.get(sequence[index : index + 3].upper(), "X")
        if residue == "*":
            break
        amino.append(residue)
    return "".join(amino)


def _gene_proteins(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (gene_id, protein) for every feature in a cluster database."""
    features = FeaturesDocument.model_validate_json(
        (path / "features.json").read_text(encoding="utf-8")
    )
    sequences = {record.id: record.sequence for record in iter_fasta(path / "sequences")}
    for locus in features.loci:
        sequence = sequences[locus.id]
        for gene in locus.genes:
            nucleotide = sequence[gene.start - 1 : gene.end]
            if gene.strand == "-":
                nucleotide = nucleotide.translate(_COMPLEMENT)[::-1]
            protein = _translate(nucleotide)
            if len(protein) >= 30:
                yield gene.gene_id, protein


def _write_fasta(path: Path, proteins: list[tuple[str, str]]) -> dict[str, str]:
    """Write numeric BLAST ids and return blast_id -> gene_id."""
    mapping: dict[str, str] = {}
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for index, (gene_id, protein) in enumerate(proteins, start=1):
            blast_id = f"g{index:07d}"
            mapping[blast_id] = gene_id
            handle.write(f">{blast_id}\n{protein}\n")
    return mapping


def _run(argv: list[str]) -> None:
    result = run_tool(argv)
    if result.returncode != 0:
        raise RuntimeError(
            f"{argv[0]} failed ({result.returncode}): {result.stderr.strip()[-800:]}"
        )


def _parse_hits(path: Path, mapping: dict[str, str]) -> list[HomologyHit]:
    hits: list[HomologyHit] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) != 7:
            continue
        query, subject = fields[0], fields[1]
        query_gene = mapping.get(query)
        subject_gene = mapping.get(subject)
        if query_gene is None or subject_gene is None:
            continue
        hits.append(
            HomologyHit(
                query=query_gene,
                subject=subject_gene,
                identity_pct=float(fields[2]),
                alignment_len=int(fields[3]),
                query_len=int(fields[4]),
                subject_len=int(fields[5]),
                evalue=float(fields[6]),
            )
        )
    return hits


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("db_dir", type=Path, help="cluster database directory")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-coverage", type=float, default=80.0)
    parser.add_argument("--evalue", type=float, default=1e-5)
    args = parser.parse_args(argv)

    db_dir: Path = args.db_dir
    output: Path = args.output or (db_dir / GENE_FAMILY_FILENAME)
    features = FeaturesDocument.model_validate_json(
        (db_dir / "features.json").read_text(encoding="utf-8")
    )
    proteins = list(_gene_proteins(db_dir))
    if not proteins:
        raise RuntimeError(f"no translatable genes found in {db_dir}")
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="gspec-families.") as workdir:
        work = Path(workdir)
        fasta = work / "genes.faa"
        mapping = _write_fasta(fasta, proteins)
        database = work / "genes.db"
        _run(["makeblastdb", "-in", str(fasta), "-dbtype", "prot", "-out", str(database)])
        hits_path = work / "hits.tsv"
        _run(
            [
                "blastp",
                "-query",
                str(fasta),
                "-db",
                str(database),
                "-outfmt",
                "6 qseqid sseqid pident length qlen slen evalue",
                "-evalue",
                str(args.evalue),
                "-max_target_seqs",
                "50",
                "-max_hsps",
                "1",
                "-num_threads",
                str(args.threads),
                "-out",
                str(hits_path),
            ]
        )
        raw_hits = _parse_hits(hits_path, mapping)
    hits = filter_homology_hits(
        raw_hits,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        max_evalue=args.evalue,
    )
    assignments = cluster_families((gene_id for gene_id, _ in proteins), hits)
    document = build_family_document(features, assignments)
    write_gene_families(document, output)
    elapsed = time.perf_counter() - started
    print(
        f"wrote {output}: {len(document.families)} families over"
        f" {len(assignments)} genes, {len(document.containment)} containment edges"
        f" in {elapsed:.1f}s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
