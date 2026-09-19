# Screening reads (FASTQ)

`gapit screen --r1/--r2` maps raw FASTQ reads against a gene database with minimap2 and calls
gene presence from alignment breadth. abricate cannot screen raw reads at all; this mode is a
gapit extension, so its defaults differ from the contig pipeline.

Requires no BLAST indexing: minimap2 indexes the db `sequences` FASTA in memory. Databases:
[./databases.md](./databases.md). Contig mode: [./screen.md](./screen.md).

## Invocation

```text
gapit screen --r1 R1[,R1b,...] [--r2 R2[,R2b,...]] --db NAME --read-type sr|map-ont|map-hifi
```

- `--r1` and `--r2` take comma-separated file lists, one entry per lane.
- Lane i pairs `r1[i]` with `r2[i]`, so the `--r2` count must equal the `--r1` count.
- All lanes aggregate into one sample: per-gene metrics pool every lane's alignments, and a
  gene can reach the breadth threshold through the union of sub-threshold lanes.
- Reads mode and positional contig files are mutually exclusive; `--r2` without `--r1` is a
  usage error (exit 2). Gzipped FASTQ works.

## Options

Reads mode runs through the same `gapit screen` command; these are the flags that apply
(transcribed from `gapit screen --help`, gapit 0.1.0). Contig-mode flags not listed here
(`--minid`, `--mincov`, `--jobs`, `--fofn`, `--noheader`, `--nopath`, `--csv`) do not apply.

| Flag | Type | Default | Description |
|---|---|---|---|
| `--r1` | str | required | Comma-separated FASTQ R1 file(s), one per lane. |
| `--r2` | str | none | Comma-separated mate FASTQ file(s); must match `--r1` count. |
| `--read-type` | sr\|map-ont\|map-hifi | `sr` | minimap2 preset for reads mode. |
| `--min-breadth` | float | `90.0` | Minimum %breadth for presence. |
| `--db` | str | `ncbi` | Database to screen against (datadir subdir). |
| `--datadir` | path | `$GAPIT_DATADIR`, then `~/.local/share/gapit/db` | Database directory. |
| `--threads` | int | `1` | minimap2 worker threads. |
| `--quiet` | flag | off | Silence stderr diagnostics. |
| `--debug` | flag | off | Verbose stderr diagnostics; echoes the minimap2 command line. |
| `--format` | tsv\|csv\|json\|md | `json` | Output format. `tsv` and `csv` are rejected in reads mode. |

`--minid` and `--mincov` do NOT apply to reads mode. Presence is decided by breadth only.

## Choosing `--read-type`

| Sequencer | Preset |
|---|---|
| Illumina (short reads, paired or single) | `sr` |
| Oxford Nanopore | `map-ont` |
| PacBio HiFi | `map-hifi` |

The default is `sr`; set the preset explicitly for long reads or mapping quality suffers.

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
you may also want to lower it. There is no per-read identity or MAPQ filter in v1.

## Output

The default is JSON, schema `gapit.reads/1`:

- `files[]` mirrors the input: each entry lists the reads it screened and the genes found.
- Gene entries sort by `breadth_pct` descending, then gene name.
- Introspect the schema with `gapit schema reads`; details in
  [./outputs.md](./outputs.md).

`--format md` produces the Markdown form with YAML frontmatter. `--format tsv` and `--format
csv` are rejected: reads results are nested per sample, not flat rows, so there is no
abricate-shaped table to emit. The refusal is a usage error, exit 2, with the usual envelope.

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
    "version": "0.1.0"
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
tool: gapit 0.1.0
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

## Notes

- **Primary alignments only.** Each read contributes at most one alignment per gene
  (`tp:A:P` records; records without a `tp` tag count too). Secondary alignments never inflate
  breadth or depth.
- **Homologous gene families.** Closely related alleles compete for the same reads. Primary-only
  assignment can misattribute shared reads between near-identical family members, so breadth for
  one member may look low when the family as a whole is well covered. There is no SNP-level
  allele calling in v1.
- **Very short genes.** With the `sr` preset, soft-clipping at alignment ends caps achievable
  breadth; genes under about 100 nt may never reach 90%. Reads mode targets normal-length genes
  (hundreds of nt and up); lower `--min-breadth` for tiny dbs.
- **Gzip input.** `.fq.gz` files work; the extension is detected the same way as contig mode.
