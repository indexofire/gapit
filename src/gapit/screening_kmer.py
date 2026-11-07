"""K-mer-only cluster screening: cluster databases accept FASTQ reads.

When a cluster database carries a ``gapit.kmerdb/1`` sidecar, reads-mode
screening bypasses minimap2 entirely — the query's canonical k-mers are
counted from FASTQ, per-locus containment provides the coverage/identity
evidence, and the differential sets deliver SNP-level ``key_kmers`` calls.
No assembly needed; the same ``gapit.typing/3`` rules evaluate identically.
"""

from __future__ import annotations

from pathlib import Path

from gapit.bundled import resolve_screen_datadir
from gapit.cluster import ClusterReport, LocusCall, load_features, load_typing
from gapit.engines import OutputFormat
from gapit.errors import usage_fail
from gapit.gbfeatures import FeaturesDocument
from gapit.kmerdb import (
    KmerDbDocument,
    containment_to_identity,
    count_fastq_kmers,
    load_kmer_set,
    locus_containment,
    read_kmerdb,
    set_hit_fraction,
)
from gapit.screening import find_database


def screen_cluster_kmer(
    r1_paths: list[Path],
    r2_paths: list[Path] | None,
    db_name: str,
    datadir: Path | None,
    min_cluster_cov: float,
    output_format: OutputFormat | None,
    quiet: bool = True,
) -> str:
    """K-mer containment screening of FASTQ reads against a cluster db
    with a kmer substrate. Produces the same TSV output as the alignment
    path (minus per-gene verdicts — no alignment happens)."""
    database = find_database(resolve_screen_datadir(datadir, db_name), db_name, quiet=quiet)
    kmer_doc = read_kmerdb(database.path)
    if kmer_doc is None:
        usage_fail(
            f"cluster database '{db_name}' has no kmer.json sidecar;"
            " k-mer reads screening requires gspec db build -K"
        )
    features = load_features(database)
    typing_document = load_typing(database)
    cols = [
        "FILE",
        "BEST_LOCUS",
        "TYPE",
        "PHENOTYPE",
        "COVERAGE",
        "IDENTITY",
        "PRESENT",
        "PARTIAL",
        "MISSING_IDS",
    ]
    lines: list[str] = ["\t".join(cols)]
    for idx, r1 in enumerate(r1_paths):
        r2 = r2_paths[idx] if r2_paths else None
        reads = [r1] + ([r2] if r2 else [])
        if not quiet:
            import click

            click.echo(f"K-mer screening: {', '.join(str(p) for p in reads)}", err=True)
        query_kmers = count_fastq_kmers(*reads, k=kmer_doc.k)
        loci = locus_calls_from_kmers(query_kmers, database.path, features, kmer_doc)
        best_locus, best_type, best_phen = "-", "-", "-"
        cov_str, id_str = "0.00", "0.00"
        if loci and loci[0].coverage_pct >= min_cluster_cov:
            top = loci[0]
            best_locus = top.locus
            best_type = top.type
            cov_str = f"{top.coverage_pct:.2f}"
            id_str = f"{top.identity_pct:.2f}"
            if typing_document is not None:
                best_phen = best_phenotype(loci, typing_document)
        lines.append(f"{r1}\t{best_locus}\t{best_type}\t{best_phen}\t{cov_str}\t{id_str}\t0\t0\t-")
    return "\n".join(lines) + "\n"


def locus_calls_from_kmers(
    query_kmers: set[int],
    db_dir: Path,
    features: FeaturesDocument,
    kmer_doc: KmerDbDocument,
) -> list[LocusCall]:
    """Per-locus LocusCall from k-mer containment: coverage = containment,
    identity = containment**(1/k) (mash ANI estimate), genes empty."""
    calls: list[tuple[float, LocusCall]] = []
    for locus in features.loci:
        cov = locus_containment(query_kmers, db_dir, locus.id)
        if cov < 0.05:
            continue
        ident = containment_to_identity(cov, kmer_doc.k)
        kmer_hits: dict[str, float] = {}
        for set_name, spec in kmer_doc.sets.items():
            if spec.of == locus.id:
                ref = load_kmer_set(db_dir, set_name)
                kmer_hits[set_name] = set_hit_fraction(query_kmers, ref)
        call = LocusCall(
            locus=locus.id,
            label=locus.id,
            type=locus.id,
            coverage_pct=round(cov * 100.0, 2),
            identity_pct=round(ident, 2),
            rank=0,
            genes=(),
            missing=(),
            kmer_hits=kmer_hits,
        )
        calls.append((cov, call))
    calls.sort(key=lambda pair: -pair[0])
    return [call.model_copy(update={"rank": i + 1}) for i, (_, call) in enumerate(calls)]


def best_phenotype(loci: list[LocusCall], typing_document: object) -> str:
    """Evaluate typing schemes on the k-mer locus calls; return the first
    scheme's phenotype (the common single-scheme case)."""
    from gapit.typing_cluster2 import evaluate_cluster2
    from gapit.typing_models import TypingDocument

    assert isinstance(typing_document, TypingDocument)
    report = ClusterReport(file="kmer", best=None, loci=tuple(loci))
    phenotypes = evaluate_cluster2(report, typing_document)
    for _, call in phenotypes.items():
        if call.phenotype is not None:
            return call.phenotype
        if call.score > 0:
            return f"ambiguous ({call.score:.3f})"
    return "-"
