# Multi-backend typing architecture

This guide explains how gapit combines alignment and k-mer backends in one
typing system, and how to design custom databases that exploit the right
backend for each discrimination problem. Read [Custom databases](./custom-db.md)
for the build mechanics; this page is the architectural companion.

## The big picture

```text
                        gapit screen
                            │
              ┌─────────────┼─────────────┐
              │             │             │
         FASTA input   FASTQ input   FASTQ input
              │             │             │
         ┌────┴────┐   ┌────┴────┐   ┌────┴────┐
         │  gene   │   │  gene   │   │ cluster │
         │   db    │   │   db    │   │ +kmer.json│
         │ BLASTn  │   │minimap2 │   │ k-mer   │
         └────┬────┘   │  (sr)   │   │containment│
              │        └────┬────┘   └────┬────┘
              │             │             │
              └─────────────┼─────────────┘
                            │
                     gapit typing
                            │
                 ┌──────────┼──────────┐
                 │          │          │
            exact_set  weighted_  cluster_match
                       genes      (cov+id+key_genes
                                   +key_kmers)
```

Two database kinds, two input kinds, three algorithm paths — every
combination works:

| Database kind | FASTA (assembly) | FASTQ (reads) |
|---|---|---|
| **gene** | BLASTn alignment | minimap2 `-x sr` alignment |
| **cluster** | minimap2 `-x asm20` alignment + k-mer hit fractions | k-mer containment (no alignment) |

The **discrimination power** varies by path:

| Discrimination problem | Alignment | K-mer | Winner |
|---|---|---|---|
| Gene presence/absence | ✅ | ✅ | Alignment (gene calls, coverage) |
| Locus family assignment | ✅ coverage/identity | ✅ containment | Either |
| **SNP-level pair** (6A vs 6B) | ❌ (~0.001 identity gap) | ✅ (1.0 vs 0.0) | **key_kmers** |
| Multi-copy gene resolution | BLASTn `--count-same-copy` | N/A | Alignment |
| Split-locus assembly | `merge_fragments` | Natural (k-mers are additive) | K-mer |

## Database kinds

### Gene databases

A flat panel of reference sequences (amplicons, genes, proteins). The
screen reports per-gene presence/absence with coverage/identity. Typing
rules are boolean:

- `exact_set`: all `requires` present + all `excludes` absent → 1.0
- `weighted_genes`: each gene carries a weight; the sum is the score

**Best for**: PCR-style marker panels (Penner multiplex, Doumith
serogrouping), virulence factor panels, small gene clusters.

**Control gene**: gene databases support reads-mode relative-depth gating
via a single-copy internal control (`lpxA`, `tuf`). See
[Reads](./reads.md).

### Cluster databases

Full genomic loci (10–50 kb) with CDS annotations. The screen reports
per-locus coverage/identity/gene calls via minimap2 alignment; typing
rules score over these:

- `cluster_match`: weighted sum of coverage + identity + `key_genes` +
  `key_kmers` components per locus

**Best for**: polysaccharide biosynthesis loci (capsule, LPS, O-antigen),
pathogenicity islands, any large variable region where gene content
matters.

**Profile**: `cluster.json` configures the coverage mode
(`primary_secondary` for multi-homologous loci) and quality gates.

## Algorithm backends in detail

### 1. BLASTn (gene, FASTA)

The workhorse: per-gene best HSP, coverage/identity thresholds, optional
`--merge-fragments` for split assemblies. Reliable, well-understood, the
abricate-compatible baseline.

### 2. minimap2 `-x sr` (gene, FASTQ)

Short-read alignment for gene panels. Depth/breadth per gene; the
`screen.json` reads section configures MAPQ, depth, and
`min_relative_depth` against a control gene for false-positive filtering.

### 3. minimap2 `-x asm20 --cs` (cluster, FASTA)

Assembly-to-locus alignment: union coverage across contigs, per-gene
verdicts (present/partial/absent), insertions/deletions map. The
`primary_secondary` coverage mode includes secondary alignments (essential
for multi-homologous loci like the Penner CPS complex).

### 4. K-mer containment (cluster, FASTQ)

No alignment at all. The query's canonical k-mers (2-bit packed, rolling
count) are intersected with each locus's full k-mer set:

- **coverage** = fraction of locus k-mers present in the query
- **identity** = `containment^(1/k)` (mash-style ANI estimate)
- **key_kmers** = fraction of the locus's *differential* set present

This path eliminates competitive mapping entirely (k-mer presence is set
membership, not best-alignment) and works on raw reads.

### 5. K-mer differential sets (cluster, both FASTA and FASTQ)

The k-mer granularity analogue of `key_genes`: the k-mers unique to one
locus against its family siblings.

```json
{"key_kmers": {"weight": 0.25, "k": 31, "set": "06A.diff"}}
```

For a SNP-level pair (6A vs 6B differ by one nucleotide), the
differential set contains the ~31 k-mers spanning the SNP:

- 6A query: `key_kmers(06A) = 1.0`, `key_kmers(06B) = 0.0`
- 6B query: `key_kmers(06B) = 1.0`, `key_kmers(06A) = 0.0`

Binary separation — no ambiguity margin can swallow this.

## Building a custom database: decision guide

| Your discrimination problem | Database kind | Typing model | Key component |
|---|---|---|---|
| PCR-style markers (presence/absence) | gene | `exact_set` | — |
| Graded gene contribution | gene | `weighted_genes` | — |
| Large locus, gene content varies | cluster | `cluster_match` | `key_genes` |
| Large locus, SNP-level pairs | cluster + kmer.json | `cluster_match` | `key_kmers` |
| Large locus, both gene content + SNPs | cluster + kmer.json | `cluster_match` | `key_genes` + `key_kmers` |
| Reads-only workflow, no assembly | cluster + kmer.json | `cluster_match` | `key_kmers` (containment path) |

### Step-by-step: cluster database with k-mer substrate

```bash
# 1. Prepare the locus FASTA (or GBK/GFF3 with CDS annotations)
# 2. Design the typing rules (gapit.typing/3)
# 3. Declare the k-mer sets (gapit.kmerdb/1)
# 4. Build

gapit db build my_sero loci.gff3 \
  --typing typing.json \
  --kmerdb kmer.json \
  --description "My organism serotyping"

# 5. Screen (FASTA)
gapit screen -d my_sero -C 60 -G 90 assembly.fna

# 6. Screen (FASTQ — no assembly needed)
gapit screen -d my_sero --r1 reads_R1.fastq --r2 reads_R2.fastq
```

The `kmer.json` declares one k size and named sets:

```json
{
  "schema": "gapit.kmerdb/1",
  "k": 31,
  "sets": {
    "typeX.diff": {"of": "typeX", "against": ["typeY", "typeZ"]}
  }
}
```

At build time gapit computes the full k-mer set of every locus and the
differential set of every declaration, storing them as sorted uint64
arrays under `kmers/`. The typing rules reference set names; the engine
validates them at build time (TYPING_UNKNOWN_KMERSET).

### Weight design for multi-path rules

Rules that must work on both FASTA and FASTQ should account for the fact
that reads-mode produces no gene calls (`key_genes` scores 0):

```json
{
  "coverage":  {"weight": 0.45, "floor": 60},
  "identity":  {"weight": 0.20, "floor": 90},
  "key_genes": {"weight": 0.10, "genes": ["geneA", "geneB"]},
  "key_kmers": {"weight": 0.25, "k": 31, "set": "typeX.diff"}
}
```

- FASTA path max: 0.45 + 0.20 + 0.10 + 0.25 = **1.00**
- FASTQ path max: 0.45 + 0.20 + 0 + 0.25 = **0.90**
- Set `cutoff` ≤ 0.75 so both paths can reach a confident call

## Real-world examples

All databases below are installed and validated; their source scripts are
archived under `~/research/`.

| Database | Organism | Kind | Rules | k-mer | SNP-level | Input |
|---|---|---|---|---|---|---|
| `cje_cps` | *C. jejuni* | gene | 38 `exact_set` | — | — | FASTA+FASTQ |
| `cje_cpsv2` | *C. jejuni* | gene | 31 `exact_set` (full loci) | — | — | FASTA+FASTQ |
| `cje_cpsv2c` | *C. jejuni* | cluster | 31 `cluster_match` | — | — | FASTA |
| `lpn_sero` | *L. pneumophila* | cluster | 13 `cluster_match` | — | — | FASTA |
| `lpn_vf` | *L. pneumophila* | gene (prot) | — (screen-only) | — | — | FASTA |
| `str_lancefield` | *Streptococcus* | gene | 5 `exact_set` + tuf | — | — | FASTA+FASTQ |
| `spn_sero` | *S. pneumoniae* | cluster | 100 `cluster_match` | ✅ k=31 | ✅ | **FASTA+FASTQ** |

### spn_sero: the full multi-backend showcase

100 cps locus references from SeroBA, 100 single-serotype rules with
differential k-mer sets. The same `gapit.typing/3` document evaluates
identically on assemblies (alignment path: coverage + identity + key_genes
+ key_kmers) and reads (k-mer containment path: estimated coverage/identity
+ key_kmers). SNP-level pairs (6A/6B, 9A/9V, 15B/15C, 18B/18C…) resolve
on both paths.

## Provenance rules (mandatory)

Every custom database must carry source metadata in its manifest and a
brief description in `info.json`. See [Custom databases → Provenance](./custom-db.md#provenance-and-license-metadata-mandatory)
for the full rule.

## Further reading

- [Custom databases](./custom-db.md) — build mechanics, provenance
- [Typing](./typing.md) — rule models in depth
- [Reads](./reads.md) — FASTQ screening for gene databases
- [Databases](./databases.md) — bundled/provider databases, naming
