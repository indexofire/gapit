# Screening reads and assemblies (FASTQ / FASTA)

`gapit screen --r1/--r2` maps raw FASTQ reads, or whole assembly FASTA files, against a gene
database with minimap2 and calls gene presence from alignment breadth. abricate cannot screen
reads at all; this mode is a gapit extension, so its defaults differ from the contig pipeline.

Requires no BLAST indexing: minimap2 indexes the db `sequences` FASTA in memory. Databases:
[./databases.md](./databases.md). Contig mode: [./screen.md](./screen.md).

## Invocation

```text
gapit screen --r1 R1[,R1b,...] [--r2 R2[,R2b,...]] --db NAME [--read-type sr|map-ont|map-hifi]
```

- `--r1` takes FASTQ reads or a FASTA assembly (see below); comma-separated, one entry per lane.
- Lane i pairs `r1[i]` with `r2[i]`, so the `--r2` count must equal the `--r1` count.
- All lanes aggregate into one sample: per-gene metrics pool every lane's alignments, and a
  gene can reach the breadth threshold through the union of sub-threshold lanes.
- Reads mode and positional contig files are mutually exclusive; `--r2` without `--r1` is a
  usage error (exit 2). Gzipped input works.
- **Input detection.** Each `--r1`/`--r2` file is detected from content at validation time:
  the first non-whitespace byte `>` means FASTA, `@` means FASTQ (gzip-wrapped files are
  peeked through the decompressor). Anything else is an input error, exit 5, code
  `INVALID_READS_FORMAT`. Mixing FASTA and FASTQ within one `--r1` list, pairing FASTA with
  `--r2`, or giving FASTA an explicit `sr`/`map-hifi` preset are usage errors (exit 2).

## Options

Reads mode runs through the same `gapit screen` command; these are the flags that apply
(transcribed from `gapit screen --help`, gapit 0.2.2). Contig-mode flags not listed here
(`--minid`, `--mincov`, `--jobs`, `--fofn`, `--noheader`, `--nopath`) do not apply.

| Flag | Type | Default | Description |
|---|---|---|---|
| `--r1` | str | required | Reads or assembly FASTA file(s), comma-separated, one per lane. |
| `--r2` | str | none | Comma-separated mate FASTQ file(s); must match `--r1` count. |
| `--read-type` | sr\|map-ont\|map-hifi | `sr` for FASTQ, `map-ont` for FASTA | minimap2 preset; resolved from the detected input when omitted. |
| `--min-breadth` | float | `90.0` | Minimum %breadth for presence. |
| `--min-identity` | float | `0.0` (off) | Minimum %identity per alignment (0–100); any nonzero value turns on [gapit.reads/2](#filtering-alignments-by-identity-and-mapq-gapitreads2) filtering. |
| `--min-mapq` | int | `0` (off) | Minimum MAPQ per alignment; any nonzero value turns on gapit.reads/2 filtering. |
| `--db` | str | `ncbi` | Database to screen against (datadir subdir). |
| `--datadir` | path | `$GAPIT_DATADIR`, then `~/.local/share/gapit/db` | Database directory. |
| `--threads` | int | `1` | minimap2 worker threads. |
| `--quiet` | flag | off | Silence stderr diagnostics (including the assembly-FASTA note). |
| `--debug` | flag | off | Verbose stderr diagnostics; echoes the minimap2 command line. |
| `--format` | tsv\|csv\|json\|md | `json` | Output format. `tsv` and `csv` are rejected in reads mode. |

`--minid` and `--mincov` do NOT apply to reads mode. Presence is decided by breadth only.

## Choosing `--read-type`

| Sequencer | Preset |
|---|---|
| Illumina (short reads, paired or single) | `sr` |
| Oxford Nanopore | `map-ont` |
| PacBio HiFi | `map-hifi` |

When `--read-type` is omitted, the preset is resolved from the detected input: `sr` for FASTQ
(the historical default), `map-ont` for assembly FASTA, announced by one stderr note
(`assembly FASTA detected; using map-ont`). An explicit preset that contradicts the input is a
usage error: assembly FASTA requires `map-ont` (exit 2). Any explicit preset is accepted for
FASTQ; set it for long reads or mapping quality suffers.

## Metrics

For each gene (the db sequence, length `tlen`), over the primary alignments of mapped reads:

| Metric | Meaning |
|---|---|
| `breadth_pct` | `100 * covered_bases / tlen`, where covered bases are the union of all primary alignments' subject spans. |
| `mean_depth` | `sum(per-base depth) / tlen`, averaged over the full gene length. |
| `reads_mapped` | Count of distinct read names with at least one primary alignment on the gene. |

A gene is `present` when `breadth_pct >= --min-breadth` (default 90.0, srst2-style). Genes with
zero mapped reads are omitted from the output entirely.

Tuning advice: the default 90% breadth is strict on purpose. Screening tiny references (under
roughly 100 nt) with `sr` rarely reaches it, because minimap2 soft-clips a few bases at each
alignment end, so calibrate `--min-breadth` downward for small custom dbs. For noisy long reads
you may also want to lower it. Per-read identity and MAPQ filtering are opt-in via
`--min-identity` / `--min-mapq` (next section).

## Output

The default is JSON, schema `gapit.reads/1`:

- `files[]` mirrors the input: each entry lists the reads it screened and the genes found.
- Gene entries sort by `breadth_pct` descending, then gene name.
- Introspect the schema with `gapit schema reads`; details in
  [./outputs.md](./outputs.md).

`--format md` produces the Markdown form with YAML frontmatter. `--format tsv` and `--format
csv` are rejected: reads results are nested per sample, not flat rows, so there is no
abricate-shaped table to emit. The refusal is a usage error, exit 2, with the usual envelope.

## Filtering alignments by identity and MAPQ (gapit.reads/2)

Reads-mode presence is breadth-only by default, which over-calls homologous gene families:
reads from a novel allele pile onto every similar db entry as low-identity primary alignments,
and breadth accumulates until the wrong gene crosses the threshold. The KP benchmark made the
failure concrete — blastn contig screening confirmed 18 genes, the reads pipeline called 11
with Illumina `sr` (split alignments starve both family members) and 140 with ONT (noisy
reads over-call nearly everything). `--min-identity` and `--min-mapq` are the opt-in fix:

- `--min-identity FLOAT` (0–100, default 0 = off): keep only alignments with
  `identity >= threshold`. Per-alignment identity is `100 * (alen - nm) / alen` over the PAF
  block length and the `NM:i:` mismatch count; a row without NM counts as 100% (cannot assess).
- `--min-mapq INT` (default 0 = off): keep only alignments with PAF MAPQ >= threshold — the
  tool for multi-mapping reads that spread breadth across repeated gene copies.
- Both filters run after the primary-only rule and before aggregation; breadth, depth and
  `reads_mapped` are computed over the surviving alignments only. Genes that lose every
  alignment drop out of the output entirely.
- With both thresholds off, nothing changes: the run stays `gapit.reads/1`, byte-identical,
  and the minimap2 invocation is untouched. Any nonzero value switches the output to the
  **`gapit.reads/2`** document — same shape, plus `params.min_identity` / `params.min_mapq`
  and a per-gene `mean_identity_pct` (alignment-length-weighted mean identity of the kept
  alignments). Introspect it with `gapit schema reads2`. Under `--format md` the frontmatter
  gains the two thresholds and each gene row gains an `Identity%` column.
- One subtlety: the `/2` invocation passes minimap2 `--cs` (the only PAF-side flag that emits
  `NM:i:`), and minimap2's emitted spans can differ slightly between the two geometries —
  expect small breadth differences between an unfiltered `/1` run and a filtered `/2` run of
  the same reads.
- The flags are reads-engine-only: combining them with the blastn contig pipeline is a usage
  error (exit 2).

Guidance: `--min-identity 95` approximates allele-level stringency (alignments from the true
gene survive at ~98-100%, homologs at ~85-92% drop out); lower it toward 90 for raw ONT reads,
whose true alignments are noisier. Use `--min-mapq 20` or higher when the db contains repeats
or duplicated gene copies and breadth splits between them.

### Worked example: the homolog fixture

The repo fixtures reproduce the failure and the fix deterministically
(`tests/data/reads2_db/homologs` + `tests/data/reads2`): the db holds a 600 nt `geneA` plus a
400 nt partial homolog `geneB` (~85% identical, absent from the sample); the sample carries
`geneA` and a novel B-like allele whose reads land on `geneB` at ~90% identity. ONT-style
reads, unfiltered — both genes are called present:

```console
$ mkdir -p /tmp/gapit-demo/readdb
$ cp -r tests/data/reads2_db/homologs /tmp/gapit-demo/readdb/
$ export GAPIT_DATADIR=/tmp/gapit-demo/readdb
$ cd tests/data/reads2
$ gapit screen --r1 ont_homologs.fq --db homologs --read-type map-ont --quiet
{
  "schema": "gapit.reads/1",
  ...
  "files": [
    {
      "reads": [
        "ont_homologs.fq"
      ],
      "genes": [
        {
          "gene": "geneB",
          ...
          "tlen": 400,
          "breadth_pct": 98.75,
          "mean_depth": 14.03,
          "reads_mapped": 15,
          "present": true
        },
        {
          "gene": "geneA",
          ...
          "tlen": 600,
          "breadth_pct": 98.5,
          "mean_depth": 14.43,
          "reads_mapped": 15,
          "present": true
        }
      ]
    }
  ]
}
```

`geneB` is a false call — the sample does not carry it. Adding `--min-identity 95` drops every
~89%-identity alignment on `geneB` (the gene vanishes; zero-read genes are omitted) and keeps
the true gene's ~98% alignments:

```console
$ gapit screen --r1 ont_homologs.fq --db homologs --read-type map-ont --min-identity 95 --quiet
{
  "schema": "gapit.reads/2",
  ...
  "params": {
    "db": "homologs",
    "read_type": "map-ont",
    "min_breadth": 90.0,
    "threads": 1,
    "min_identity": 95.0,
    "min_mapq": 0
  },
  "files": [
    {
      "reads": [
        "ont_homologs.fq"
      ],
      "genes": [
        {
          "gene": "geneA",
          ...
          "breadth_pct": 100.0,
          "mean_depth": 14.98,
          "reads_mapped": 15,
          "present": true,
          "mean_identity_pct": 98.04
        }
      ]
    }
  ]
}
```

The short-read fixture behaves the same with one calibration: `sr` soft-clipping caps
`geneB`'s unfiltered breadth at ~89.8%, so the unfiltered leg runs with `--min-breadth 80`
(`geneA` 94.83% / `geneB` 89.75%, both present; filtered: `geneA` 96.5% at identity 100.0,
`geneB` gone). The Markdown form of a filtered run:

```console
$ gapit screen --r1 ont_homologs.fq --db homologs --read-type map-ont --min-identity 95 --format md --quiet
---
schema: gapit.reads/2
tool: gapit 0.2.2
created_at: 2026-09-20T14:23:33Z
db: homologs
read_type: map-ont
min_breadth: 90.0
threads: 1
min_identity: 95.0
min_mapq: 0
files: 1
genes_found: 1
---

# gapit read screening report

## `ont_homologs.fq`

| Gene | Breadth% | Depth | Reads | Present | Database | Accession | Product | Resistance | Identity% |
|---|---|---|---|---|---|---|---|---|---|
| geneA | 100.00 | 14.98 | 15 | yes | homologs | SYN-A | true allele carried by the sample | TETRACYCLINE | 98.04 |
```

## Worked example

The repo test fixtures include a tiny reads db (`tinyreads`, one 522 nt `tetX` gene plus a
partial-coverage `sulY`) and small FASTQ files. Copy the db to a temp datadir and run from the
reads directory:

```console
$ mkdir -p /tmp/gapit-demo/readdb
$ cp -r tests/data/reads_db/tinyreads /tmp/gapit-demo/readdb/
$ export GAPIT_DATADIR=/tmp/gapit-demo/readdb
$ cd tests/data/reads
$ gapit screen --r1 tetx_full.fq --db tinyreads
Screening reads: tetx_full.fq
Detected 1 present genes in tetx_full.fq
{
  "schema": "gapit.reads/1",
  "tool": {
    "name": "gapit",
    "version": "0.2.2"
  },
  "created_at": "2026-09-19T01:12:25Z",
  "params": {
    "db": "tinyreads",
    "read_type": "sr",
    "min_breadth": 90.0,
    "threads": 1
  },
  "files": [
    {
      "reads": [
        "tetx_full.fq"
      ],
      "genes": [
        {
          "gene": "tetX",
          "database": "tinyreads",
          "accession": "SYN-001",
          "product": "extended resistance determinant tetX",
          "resistance": "TETRACYCLINE",
          "tlen": 522,
          "breadth_pct": 97.7,
          "mean_depth": 2.09,
          "reads_mapped": 12,
          "present": true
        }
      ]
    }
  ]
}
```

The `Screening reads:` and `Detected N present genes` lines are stderr; stdout is pure JSON.

### Paired-end and multiple lanes

One paired lane: pass both mates. Two single-end lanes as one sample: comma-join them. Both
aggregate to the same result here:

```console
$ gapit screen --r1 tetx_R1.fq --r2 tetx_R2.fq --db tinyreads --quiet
{
  "schema": "gapit.reads/1",
  ...
  "files": [
    {
      "reads": [
        "tetx_R1.fq",
        "tetx_R2.fq"
      ],
      "genes": [
        {
          "gene": "tetX",
          ...
          "breadth_pct": 97.7,
          "mean_depth": 2.09,
          "reads_mapped": 12,
          "present": true
        }
      ]
    }
  ]
}
$ gapit screen --r1 tetx_lane1.fq,tetx_lane2.fq --db tinyreads --quiet
{
  "schema": "gapit.reads/1",
  ...
  "files": [
    {
      "reads": [
        "tetx_lane1.fq",
        "tetx_lane2.fq"
      ],
      "genes": [
        {
          "gene": "tetX",
          ...
          "breadth_pct": 97.7,
          "mean_depth": 2.09,
          "reads_mapped": 12,
          "present": true
        }
      ]
    }
  ]
}
```

(`...` marks output identical to the first JSON block above.)

### Markdown output

```console
$ gapit screen --r1 tetx_full.fq --db tinyreads --format md
---
schema: gapit.reads/1
tool: gapit 0.2.2
created_at: 2026-09-19T01:12:25Z
db: tinyreads
read_type: sr
min_breadth: 90.0
threads: 1
files: 1
genes_found: 1
---

# gapit read screening report

## `tetx_full.fq`

| Gene | Breadth% | Depth | Reads | Present | Database | Accession | Product | Resistance |
|---|---|---|---|---|---|---|---|---|
| tetX | 97.70 | 2.09 | 12 | yes | tinyreads | SYN-001 | extended resistance determinant tetX | TETRACYCLINE |
```

### Tuning `--min-breadth`

`suly_partial.fq` covers 63.98% of the 261 nt `sulY` gene. At the default threshold the gene is
reported with `present: false`; lowering the threshold flips the call without touching the
metrics:

```console
$ gapit screen --r1 suly_partial.fq --db tinyreads --quiet
{
  "schema": "gapit.reads/1",
  ...
  "genes": [
    {
      "gene": "sulY",
      ...
      "tlen": 261,
      "breadth_pct": 63.98,
      "mean_depth": 1.0,
      "reads_mapped": 3,
      "present": false
    }
  ]
}
$ gapit screen --r1 suly_partial.fq --db tinyreads --min-breadth 50 --quiet
{
  "schema": "gapit.reads/1",
  ...
  "genes": [
    {
      "gene": "sulY",
      ...
      "tlen": 261,
      "breadth_pct": 63.98,
      "mean_depth": 1.0,
      "reads_mapped": 3,
      "present": true
    }
  ]
}
```

### TSV and CSV are rejected

```console
$ gapit screen --r1 tetx_full.fq --db tinyreads --format tsv; echo "exit=$?"
{"schema":"gapit.error/1","code":"USAGE_ERROR","message":"--format tsv|csv is not available in reads mode (use json or md)","context":{}}
exit=2
```

## Screening assemblies (fast presence survey)

The primary spelling is `gapit screen --aligner minimap2 assembly.fa`: the positional file is
routed through the minimap2 engine, and every input must be FASTA content, gzipped or plain —
a FASTQ file there is a usage error (exit 2). The `--r1 assembly.fa` spelling is equivalent.
minimap2 takes FASTA queries natively, so gapit passes the file through exactly as it does for
FASTQ; what changes is the preset: content detection forces `map-ont` (a contiguous 522 nt
contig, for example, aligns to only ~14% of its gene under `sr` because short-read soft-clipping
wrecks long-query alignments), and gapit says so on stderr. Each contig acts as one long read:
`reads_mapped` counts contigs, `mean_depth` hovers around the covered fraction, and `present`
still means `breadth_pct >= --min-breadth`. The output stays `gapit.reads/1`;
`params.read_type` reports the resolved preset.

Using the tinyreads fixture as a stand-in assembly (any multi-contig FASTA behaves the same):

```console
$ mkdir -p /tmp/gapit-demo/readdb
$ cp -r tests/data/reads_db/tinyreads /tmp/gapit-demo/readdb/
$ export GAPIT_DATADIR=/tmp/gapit-demo/readdb
$ cp tests/data/reads_db/tinyreads/sequences /tmp/gapit-demo/assembly.fa
$ gapit screen --aligner minimap2 /tmp/gapit-demo/assembly.fa --db tinyreads
assembly FASTA detected; using map-ont
Screening reads: /tmp/gapit-demo/assembly.fa
Detected 2 present genes in /tmp/gapit-demo/assembly.fa
{
  "schema": "gapit.reads/1",
  "tool": {
    "name": "gapit",
    "version": "0.2.2"
  },
  "created_at": "2026-09-19T14:23:00Z",
  "params": {
    "db": "tinyreads",
    "read_type": "map-ont",
    "min_breadth": 90.0,
    "threads": 1
  },
  "files": [
    {
      "reads": [
        "/tmp/gapit-demo/assembly.fa"
      ],
      "genes": [
        {
          "gene": "tetX",
          "database": "tinyreads",
          "accession": "SYN-001",
          "product": "extended resistance determinant tetX",
          "resistance": "TETRACYCLINE",
          "tlen": 522,
          "breadth_pct": 97.89,
          "mean_depth": 0.98,
          "reads_mapped": 1,
          "present": true
        },
        {
          "gene": "sulY",
          "database": "tinyreads",
          "accession": "SYN-002",
          "product": "partial coverage test determinant sulY",
          "resistance": "SULFONAMIDE",
          "tlen": 261,
          "breadth_pct": 95.79,
          "mean_depth": 0.96,
          "reads_mapped": 1,
          "present": true
        }
      ]
    }
  ]
}
```

(The `assembly FASTA detected`, `Screening reads:`, and `Detected ...` lines are stderr;
`--quiet` silences the note and the chatter, leaving stdout byte-identical.)

### The two-stage pattern: survey, then confirm

Because the minimap2 stage skips BLAST indexing entirely, screening an assembly through the
minimap2 engine is roughly an order of magnitude faster than the BLAST contig pipeline, at the
cost of allele-level precision. That trade suggests a two-stage workflow over many samples:

1. **Survey** every sample with minimap2 and a deliberately relaxed breadth floor.
2. **Confirm** only the positives (or only the samples with any hit) with the regular contig
   pipeline, which applies the identity and coverage floors at abricate parity.

Real numbers, one K. pneumoniae RefSeq assembly (GCF_000240185.1, 5.3 Mb, `--db ncbi`,
single-threaded, gapit 0.2.2):

```console
$ # Stage 1: survey, ~0.9 s
$ gapit screen --aligner minimap2 kpneu_mgh78578.fna.gz --db ncbi --min-breadth 50 --quiet
{
  "schema": "gapit.reads/1",
  ...
  "params": {
    "db": "ncbi",
    "read_type": "map-ont",
    "min_breadth": 50.0,
    "threads": 1
  },
  "files": [
    {
      "reads": ["kpneu_mgh78578.fna.gz"],
      "genes": [
        { "gene": "aph(3'')-Ib", "breadth_pct": 99.75, "present": true, ... },
        { "gene": "dfrA12", "breadth_pct": 99.6, "present": true, ... },
        { "gene": "tet(G)", "breadth_pct": 99.57, "present": true, ... },
        { "gene": "blaKPC-2", "breadth_pct": 99.55, "present": true, ... },
        ... 11 more present genes (floR2 dfrA50 blaCTX-M-14 rmtB1 aadA2 sul2 blaTEM-1 aph(6)-Id fosA6 blaSHV-155 aac(3)-IId) ...
        { "gene": "sul1", "breadth_pct": 62.74, "present": true, ... },
        { "gene": "tmexD3", "breadth_pct": 60.22, "present": true, ... },
        { "gene": "tmexD2", "breadth_pct": 56.2, "present": true, ... }
      ]
    }
  ]
}
$ # Stage 2: confirm the positives with the contig pipeline, ~13.6 s
$ gapit screen kpneu_mgh78578.fna.gz --db ncbi --quiet
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
kpneu_mgh78578.fna.gz	NC_016838.1	31843	32718	-	blaCTX-M-14	1-876/876	===============	0/0	100.00	100.00	ncbi	NG_048929.1	extended-spectrum class A beta-lactamase CTX-M-14	CEPHALOSPORIN
kpneu_mgh78578.fna.gz	NC_016838.1	108291	108766	-	dfrA50	2-477/477	===============	0/0	99.79	98.74	ncbi	NG_242637.1	trimethoprim-resistant dihydrofolate reductase DfrA50	TRIMETHOPRIM
... 16 more rows: aac(3)-IId aadA2 aph(3'')-Ib aph(6)-Id blaKPC-2 blaSHV-158 blaTEM-1 dfrA12 floR2 fosA6 rmtB1 sul2 tet(G) ...
```

(`...` elides fields and rows; the survey JSON is trimmed to the genes discussed below.)

On this genome the survey ran in 0.9 s versus 13.6 s for BLAST (about 17x), and 14 of 18
confirmed genes appear in both lists. The differences are the point:

- **Family-level, not allele-level, resolution.** The survey flagged the blaSHV locus as
  `blaSHV-155`; BLAST confirms the locus but calls the allele `blaSHV-158`. minimap2's
  primary-only assignment hands a contig shared between near-identical family members to one
  of them, so treat survey gene names as family-level hints and let the confirm stage name
  alleles.
- **No identity floor.** Reads-mode presence is breadth-only. `sul1` (62.7% breadth) and the
  `tmexD2`/`tmexD3` pair (56-60%) pass the relaxed 50% survey floor but are partial or
  divergent loci that the confirm stage's 80/80 identity/coverage thresholds reject. A
  relaxed survey floor trades precision for recall on purpose; the exact filter belongs to
  stage 2.

When you need exact alleles in one pass, use the contig pipeline directly; use the survey when
you need gene-family answers from many assemblies quickly.

## Notes

- **Primary alignments only.** Each read contributes at most one alignment per gene
  (`tp:A:P` records; records without a `tp` tag count too). Secondary alignments never inflate
  breadth or depth.
- **Homologous gene families.** Closely related alleles compete for the same reads. Primary-only
  assignment can misattribute shared reads between near-identical family members, so breadth for
  one member may look low when the family as a whole is well covered, and reads from a novel
  allele can push an absent homolog over the breadth floor. `--min-identity` (previous section)
  removes the low-identity side of that failure. There is no SNP-level allele calling in v1.
- **Very short genes.** With the `sr` preset, soft-clipping at alignment ends caps achievable
  breadth; genes under about 100 nt may never reach 90%. Reads mode targets normal-length genes
  (hundreds of nt and up); lower `--min-breadth` for tiny dbs.
- **Gzip input.** `.gz` files work for FASTQ and FASTA alike; input detection peeks through the
  gzip wrapper, and minimap2 decompresses natively.
