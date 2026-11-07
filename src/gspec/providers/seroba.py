"""SeroBA S. pneumoniae cps locus provider (cluster kind, download-on-fetch).

Downloads SeroBA reference.fasta (GPL-3.0, sequences INSDC public domain),
annotates CDS with pyrodigal into a GFF3, and builds a cluster database
with gapit.typing/3 key_kmers differential sets (k=31) resolving SNP-level
serotype pairs. Both FASTA (alignment) and FASTQ (k-mer containment) paths
are supported.

GPL-3.0: nothing is bundled with gspec (MIT); `gspec db fetch spn_sero`
downloads and builds locally. Cite SeroBA when using results.
"""

import json
from pathlib import Path

from gspec.db_info_models import DbGuidance, DbReference
from gspec.providers.cluster_common import ClusterProvider

NAME = "spn_sero"
_RAW = "https://raw.githubusercontent.com/sanger-pathogens/seroba/master/database"

# Committed wheel assets (typing + kmer.json are the curator's own design;
# the sequences themselves are INSDC public domain, but the SeroBA
# compilation is GPL-3.0 so the download stays fetch-on-demand)
_TYPING_PATH = Path(__file__).parent / "data" / "spn_sero_typing.json"
_KMERDB_PATH = Path(__file__).parent / "data" / "spn_sero_kmer.json"
_DROP_NT = ("alternative_aliB_NT", "Swiss_NT")


def preprocess(workdir: Path) -> Path:
    """Annotate the downloaded FASTA with pyrodigal and write a GFF3.

    The SeroBA reference is a plain multi-FASTA (one cps locus per
    serotype). The cluster pipeline needs a GFF3 with CDS features; this
    hook runs pyrodigal on every non-NT reference and writes the result.
    """
    import pyrodigal

    fasta_path = workdir / "reference.fasta"
    gff3_path = workdir / "loci.gff3"
    orf = pyrodigal.GeneFinder(meta=True)
    gff_rows: list[str] = []
    fasta_lines: list[str] = []
    name: str | None = None
    seq_buf: list[str] = []
    sequences: dict[str, str] = {}
    for line in fasta_path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name and name not in _DROP_NT:
                sequences[name] = "".join(seq_buf)
            name = line[1:].strip()
            seq_buf = []
        else:
            seq_buf.append(line.strip())
    if name and name not in _DROP_NT:
        sequences[name] = "".join(seq_buf)
    for locus_id, seq in sequences.items():
        genes = orf.find_genes(seq.encode())
        for i, gene in enumerate(genes, start=1):
            gid = f"{locus_id}__{i:03d}"
            strand = "+" if gene.strand == 1 else "-"
            gff_rows.append(
                f"{locus_id}\tpyrodigal\tCDS\t{gene.begin + 1}\t{gene.end}\t."
                f"\t{strand}\t0\tID={gid};Name={gid}"
            )
        fasta_lines.append(f">{locus_id}")
        for j in range(0, len(seq), 60):
            fasta_lines.append(seq[j : j + 60])
    gff3_path.write_text(
        "##gff-version 3\n" + "\n".join(gff_rows) + "\n##FASTA\n" + "\n".join(fasta_lines) + "\n",
        encoding="utf-8",
    )
    # Replace the FASTA with the GFF3 so _merged_loci (which iterates the
    # source_urls basenames) picks up only the annotated GFF3.
    fasta_path.unlink()
    return gff3_path


def finalize(workdir: Path, build_dir: Path) -> None:
    """Copy kmer.json, materialize .u64 sets, and write cluster.json."""
    from gspec.kmerdb import KmerDbDocument, materialize_kmer_sets

    if _KMERDB_PATH.is_file():
        kmer_doc = KmerDbDocument.model_validate_json(_KMERDB_PATH.read_text(encoding="utf-8"))
        (build_dir / "kmer.json").write_text(
            kmer_doc.model_dump_json(indent=1, by_alias=True), encoding="utf-8"
        )
        materialize_kmer_sets(build_dir, kmer_doc, build_dir / "sequences")
    cluster_profile = {
        "schema": "gapit.cluster_profile/1",
        "coverage_mode": "primary_secondary",
        "fallback_evidence": True,
        "min_gene_cov": 50.0,
        "locus_typeable_coverage": 50.0,
        "locus_typeable_identity": 90.0,
        "locus_quality_min_coverage": 80.0,
        "locus_quality_min_identity": 90.0,
        "locus_quality_max_missing_genes": 3,
        "locus_quality_max_extra_panel_genes": 2,
    }
    (build_dir / "cluster.json").write_text(
        json.dumps(cluster_profile, indent=1) + "\n", encoding="utf-8"
    )


PROVIDER = ClusterProvider(
    name=NAME,
    description="S. pneumoniae cps serotyping (SeroBA, SNP-level k-mer, FASTA+FASTQ)",
    vendor="SeroBA (Sanger Pathogens)",
    source_urls=(f"{_RAW}/reference.fasta",),
    license="GPL-3.0 (SeroBA database; sequences INSDC public domain)",
    note=(
        "SeroBA reference.fasta (102 cps loci; NT refs excluded), pyrodigal"
        " CDS annotation, gapit.typing/3 with key_kmers (k=31) resolving"
        " SNP-level serotype pairs; FASTA and FASTQ dual path"
    ),
    typing=_TYPING_PATH if _TYPING_PATH.is_file() else None,
    preprocess=preprocess,
    finalize=finalize,
    guidance=DbGuidance(
        references=(
            DbReference(
                citation=(
                    "Epping L, van Tonder AJ, Bray JE, et al. SeroBA: rapid"
                    " high-throughput serotyping of Streptococcus pneumoniae from"
                    " whole genome sequencing data. Microb Genom 2018;4(9):e000186."
                ),
                doi="10.1099/mgen.0.000186",
                url="https://github.com/sanger-pathogens/seroba",
            ),
            DbReference(
                citation=(
                    "Bentley SD, Aanensen DM, Mavroidi A, et al. Genetic analysis"
                    " of the capsular biosynthetic locus from all 90 pneumococcal"
                    " serotypes. PLoS Genet 2006;2(3):e31."
                ),
                doi="10.1371/journal.pgen.0020031",
            ),
        ),
        caveats=(
            "key_kmers resolve SNP-level pairs (6A/6B, 9A/9V, 15B/15C);"
            " phase-variable variants may flip with switch state",
            "Rough/unencapsulated strains (R6-like) correctly report non-typeable",
            "GPL-3.0: fetched on demand, never bundled into the gspec wheel",
            "Both assembly (FASTA) and reads (FASTQ) supported",
        ),
    ),
)
