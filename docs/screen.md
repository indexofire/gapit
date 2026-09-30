# Screening contigs

`gapit screen` aligns contig files against a gene database with BLAST and reports which genes
are present, using the same pipeline and hit rules as abricate. Output is abricate-compatible
TSV by default; JSON and Markdown are first-class alternatives (see [./outputs.md](./outputs.md)).

Database setup is covered in [./databases.md](./databases.md); a full walk-through lives in
[./quickstart.md](./quickstart.md).

## Options

Transcribed from `gapit screen --help` (gapit 0.4.0). Flags marked *reads mode* apply only when
you pass `--r1`/`--r2`; they are documented in [./reads.md](./reads.md).

| Flag | Type | Default | Description |
|---|---|---|---|
| `FILE...` | path(s) | required* | Input FASTA/GBK/EMBL contig file(s) to screen. |
| `--db` | str | `ncbi` | Database to screen against (datadir subdir). |
| `--datadir` | path | `$GAPIT_DATADIR`, then `~/.local/share/gapit/db` | Database directory. |
| `--minid` | float | `80.0` | Minimum %identity, `0 < x <= 100`. Enforced inside BLAST via `-perc_identity`. |
| `--mincov` | float | `80.0` | Minimum %coverage, `0 <= x <= 100`. Post-filter on the unrounded float. |
| `--threads` | int | `1` | BLAST worker threads (passed to `-num_threads`). |
| `--jobs` | int | `1` | Screen N input files concurrently. Output order is always input order. |
| `--merge-fragments` | flag | off | Merge gene fragments split across contigs (gapit extension; see below). |
| `--fofn` | path | none | File of filenames; replaces the positional FILEs. |
| `--quiet` | flag | off | Silence stderr diagnostics. |
| `--noheader` | flag | off | Suppress the `#FILE ...` header row. |
| `--nopath` | flag | off | Basename the FILE column. |
| `--debug` | flag | off | Verbose stderr diagnostics; echoes each external command line. |
| `--format` | tsv\|csv\|json\|md | `tsv` | Output format (reads mode defaults to json). |
| `--aligner` | blastn\|minimap2 | input-based | Alignment engine (default: blastn for contig files, minimap2 for `--r1`/`--r2` reads). `--aligner minimap2` routes positional FASTA assemblies through the minimap2 engine (FASTA content required; see [./reads.md](./reads.md)). |
| `--r1` | str | none | *Reads mode.* Reads or assembly FASTA file(s), comma-separated, one per lane. |
| `--r2` | str | none | *Reads mode.* Comma-separated mate FASTQ file(s); must match `--r1` count. |
| `--read-type` | sr\|map-ont\|map-hifi | `sr` for FASTQ, `map-ont` for FASTA | *Reads mode.* minimap2 preset; resolved from the detected input when omitted. |
| `--min-breadth` | float | `90.0` | *Reads mode.* Minimum %breadth for presence. |
| `--min-gene-cov` | float | `90.0` | *Cluster databases only.* Minimum %coverage for a gene `present` verdict. |
| `--min-gene-id` | float | `90.0` | *Cluster databases only.* Minimum %identity for a gene `present` verdict. |
| `--min-cluster-cov` | float | `96.0` | *Cluster databases only.* Minimum locus %coverage for a best-locus call. |

\* Positional FILEs or `--fofn`, or reads mode via `--r1`. Positional files and `--r1`/`--r2`
are mutually exclusive.

## Input files

Normalization is native (no external `any2fasta`): plain FASTA, gzipped and bzip2-compressed
FASTA, FASTQ, GBK, and EMBL all work. The converted FASTA is buffered in memory (genome-scale
assemblies are a few MB) and fed to blastn on stdin. If normalization fails (not a sequence
file), gapit prints a `gapit.error/1` envelope on stderr and exits 5.

A `--fofn` file lists one path per line and replaces positional arguments entirely.

## Worked example

The `tinyamr` fixture db ships in the repo test data (three short AMR genes). Copy it to a
temp datadir, index it, and point `$GAPIT_DATADIR` at it. Never write into
`~/.local/share/gapit/db` for experiments.

```console
$ mkdir -p /tmp/gapit-demo/datadir/tinyamr
$ cp tests/data/db/tinyamr/sequences /tmp/gapit-demo/datadir/tinyamr/sequences
$ gapit setupdb --datadir /tmp/gapit-demo/datadir
Indexed tinyamr (3 sequences, nucl)
$ export GAPIT_DATADIR=/tmp/gapit-demo/datadir
```

### Default TSV

```console
$ gapit screen tests/data/contigs/full.fa --db tinyamr
Processing: tests/data/contigs/full.fa
Found 1 genes in tests/data/contigs/full.fa
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
tests/data/contigs/full.fa	contig1	1	79	+	tetA	1-79/79	===============	0/0	100.00	100.00	tinyamr	NC_000913.3:100-900	tetracycline efflux pump TetA	TETRACYCLINE
```

`Processing:` and `Found N genes` lines go to stderr; the header and hit rows are stdout. With
stderr discarded the output is pure data:

```console
$ gapit screen tests/data/contigs/full.fa --db tinyamr 2>/dev/null
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
tests/data/contigs/full.fa	contig1	1	79	+	tetA	1-79/79	===============	0/0	100.00	100.00	tinyamr	NC_000913.3:100-900	tetracycline efflux pump TetA	TETRACYCLINE
```

### JSON

```console
$ gapit screen tests/data/contigs/full.fa --db tinyamr --format json
{
  "schema": "gapit.report/1",
  "tool": {
    "name": "gapit",
    "version": "0.4.0"
  },
  "created_at": "2026-09-19T01:12:04Z",
  "params": {
    "db": "tinyamr",
    "minid": 80.0,
    "mincov": 80.0,
    "threads": 1
  },
  "files": [
    {
      "file": "tests/data/contigs/full.fa",
      "hits": [
        {
          "sequence": "contig1",
          "start": 1,
          "end": 79,
          "strand": "+",
          "gene": "tetA",
          "coverage": "1-79/79",
          "coverage_map": "===============",
          "gaps": "0/0",
          "coverage_pct": 100.0,
          "identity_pct": 100.0,
          "database": "tinyamr",
          "accession": "NC_000913.3:100-900",
          "product": "tetracycline efflux pump TetA",
          "resistance": "TETRACYCLINE"
        }
      ]
    }
  ]
}
```

Field names are snake_case with units explicit (`identity_pct`, `coverage_pct`). The schema is
introspectable: `gapit schema report`. See [./outputs.md](./outputs.md).

### How `--mincov` filters

`tests/data/contigs/partial.fa` carries the first 44 bases of the 88 bp `blaTEM-1` reference,
so coverage is 50%. The default threshold drops it:

```console
$ gapit screen tests/data/contigs/partial.fa --db tinyamr --quiet
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
```

Lowering the threshold keeps it:

```console
$ gapit screen tests/data/contigs/partial.fa --db tinyamr --mincov 50 --quiet
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
tests/data/contigs/partial.fa	contig1	1	44	+	blaTEM-1	1-44/88	========.......	0/0	50.00	100.00	tinyamr	J01749.1:1-861	class A beta-lactamase TEM-1	BETA-LACTAM
```

### Compressed input, CSV, path and header control

```console
$ gzip -c tests/data/contigs/gap.fa > /tmp/gapit-demo/gap.fa.gz
$ gapit screen /tmp/gapit-demo/gap.fa.gz --db tinyamr --quiet --nopath
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
gap.fa.gz	contig1	1	97	+	sul1	1-94/94	========/======	1/3	100.00	96.91	tinyamr	U12338.4:1-940	sulfonamide-resistant dihydropteroate synthase Sul1	SULFONAMIDE
$ gapit screen tests/data/contigs/full.fa --db tinyamr --format csv --quiet --nopath
#FILE,SEQUENCE,START,END,STRAND,GENE,COVERAGE,COVERAGE_MAP,GAPS,%COVERAGE,%IDENTITY,DATABASE,ACCESSION,PRODUCT,RESISTANCE
full.fa,contig1,1,79,+,tetA,1-79/79,===============,0/0,100.00,100.00,tinyamr,NC_000913.3:100-900,tetracycline efflux pump TetA,TETRACYCLINE
```

### Many files: `--fofn` and `--jobs`

```console
$ printf '%s\n' tests/data/contigs/full.fa tests/data/contigs/gap.fa > /tmp/gapit-demo/files.txt
$ gapit screen --fofn /tmp/gapit-demo/files.txt --db tinyamr --quiet
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
tests/data/contigs/full.fa	contig1	1	79	+	tetA	1-79/79	===============	0/0	100.00	100.00	tinyamr	NC_000913.3:100-900	tetracycline efflux pump TetA	TETRACYCLINE
tests/data/contigs/gap.fa	contig1	1	97	+	sul1	1-94/94	========/======	1/3	100.00	96.91	tinyamr	U12338.4:1-940	sulfonamide-resistant dihydropteroate synthase Sul1	SULFONAMIDE
$ gapit screen tests/data/contigs/full.fa tests/data/contigs/gap.fa --db tinyamr --jobs 2 --quiet
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
tests/data/contigs/full.fa	contig1	1	79	+	tetA	1-79/79	===============	0/0	100.00	100.00	tinyamr	NC_000913.3:100-900	tetracycline efflux pump TetA	TETRACYCLINE
tests/data/contigs/gap.fa	contig1	1	97	+	sul1	1-94/94	========/======	1/3	100.00	96.91	tinyamr	U12338.4:1-940	sulfonamide-resistant dihydropteroate synthase Sul1	SULFONAMIDE
```

With `--jobs 2` the two files are screened concurrently but stdout stays in input order, so the
output is identical to the sequential run.

### `--debug`

Echoes the native normalization step and every external command line to stderr:

```console
$ gapit screen tests/data/contigs/full.fa --db tinyamr --debug 2>&1 >/dev/null
Processing: tests/data/contigs/full.fa
gapit: normalize: tests/data/contigs/full.fa (fasta)
gapit: run: blastn -task blastn -dust no -perc_identity 80.0 -db /tmp/gapit-demo/datadir/tinyamr/sequences -outfmt '6 qseqid qstart qend qlen sseqid sstart send slen sstrand evalue length pident gaps gapopen stitle' -num_threads 1 -evalue 1E-20 -culling_limit 1 -max_target_seqs 10000
Found 1 genes in tests/data/contigs/full.fa
```

stdout is unchanged by `--debug`; it is a stderr-only flag.

## Fragment merging (`--merge-fragments`)

Draft assemblies often break one gene across a contig boundary: each contig then carries only
a partial fragment of the gene, every fragment fails `--mincov 80` on its own, and the gene
goes unreported. `--merge-fragments` (a gapit extension, off by default; blastn contig mode
only — it is a usage error in reads mode or with `--aligner minimap2`) rescues exactly those
genes:

- After the usual per-hit filtering, gapit groups the sub-`mincov` fragments of a gene that
  were deduped on `(contig, start, end)` and passed `--minid`.
- If the **union** of the fragments' subject intervals covers `>= mincov` of the gene
  (`100 * union_length / slen`, compared on the unrounded float, same boundary rule as the
  per-hit filter) **and** at least two fragments are involved, one merged hit is reported.
- Genes that already have one individually passing hit are never merged — their stray
  partials are ignored. Groups whose union stays below `mincov` still report nothing.

Merged-row semantics:

- `%COVERAGE` is the union coverage; `%IDENTITY` is the aligned-length-weighted mean of the
  fragment identities; `GAPS` and the `COVERAGE_MAP` '/' marker come from summed gap counts,
  with the map binned over the interval union (same 15-char arithmetic as normal rows).
- `SEQUENCE` is the comma-joined list of contig names carrying the fragments (sorted; fine
  for TSV — avoid `--format csv` for merged rows, where a comma inside a field is ambiguous).
- `START`/`END`/`STRAND` come from the **anchor**: the fragment with the largest aligned
  length (ties: first in `(contig, start)` order). `COVERAGE` is the union's bounding span.
- JSON (`gapit.report/1`) additionally carries `"merged": true` and a `fragments` array with
  each fragment's `{contig, start, end, strand, identity_pct, coverage_pct}` — additive
  optional fields, absent from every non-merged row. Markdown adds one detail line per
  merged gene under the table.

Without the flag nothing about this path runs: default output stays byte-identical to
abricate (the Notes below still apply in full).

## Cluster databases (`kind: cluster`)

When `--db` names a cluster database — built from GBK/GFF via
[`gapit db build`](./databases.md#cluster-databases-gbkgff) or fetched from a kaptive
provider — `gapit screen` dispatches to the **cluster engine** instead of the blastn gene
pipeline. Nothing about the gene path changes: cluster-kind dispatch is additive and
gene-kind screening is byte-identical to the pre-cluster behavior.

How it works:

- One **minimap2 `asm20`** invocation per input file (query = the file's contigs, target =
  the database's locus FASTA, `--cs` short-form alignment trace).
- Per-**locus** union coverage and identity come from the cs walk across all primary
  alignments — records from multiple contigs **union**, so a locus fragmented across contig
  boundaries still screens as fully present.
- Every annotated **gene** gets a verdict: `present` (≥ `--min-gene-cov` coverage AND ≥
  `--min-gene-id` identity), `partial` (≥ 50% coverage but failing a present threshold), or
  `absent`.
- Loci rank by coverage, identity, covered bases, then id; the rank-1 locus becomes the
  file's **best call** when its coverage reaches `--min-cluster-cov` (default 96, a
  kaptive-style confidence floor). Below the floor no call is made.

Guards (usage errors, exit 2): `--minid`/`--mincov`/`--merge-fragments`/`--jobs`/`--aligner`
are gene-engine options and are rejected on cluster databases; the three `--min-gene-*`/
`--min-cluster-cov` flags are cluster-only and rejected elsewhere, including reads mode
(cluster screening is assembly-FASTA only in v1).

### Phenotype calls (`typing.json`)

A cluster database may carry a `gapit.typing/1` scoring spec (installed with
`gapit db build --typing FILE`). When it does, each file's best call is annotated with a
**phenotype**:

- Rules score independently — `weighted_genes` (gene presence with an identity floor and
  optional negative markers), `cluster_match` (weighted locus coverage/identity/key-genes
  components with per-component floors), or `learned_linear` (a trained sigmoid model over
  named features) — and the decision layer applies the document's `cutoff` and
  `ambiguity_margin`: a clear winner is called with `high` confidence, two rules inside the
  margin yield an `ambiguous` call (phenotype null, both candidates listed), and a
  sub-cutoff best falls back to the document's fallback string.
- JSON (`gapit.cluster/1`) carries `best.phenotype` plus an additive `best.phenotype_detail`
  block: `{score, confidence, components[], runner_up, ambiguous[]}` — the explainable
  per-rule/per-component breakdown.
- The TSV/CSV header gains a **PHENOTYPE** column after TYPE (`-` when ambiguous or
  uncalled); untyped cluster databases keep the untyped header byte-identically.
- A rule referencing a gene/locus the database lacks fails up front with the typed
  `TYPING_UNKNOWN_GENE` error (also enforced at `db build --typing` time).

Untyped cluster databases (the kaptive fetches, for now) report locus calls only — the
kaptive-style output — with `phenotype` null.

## Notes

- **Coverage filter on the unrounded float.** `%COVERAGE = 100 * (length - gaps) / slen` is
  compared to `--mincov` before rounding, while display is `%.2f`. A hit at 79.996% prints as
  `80.00` but fails the default `--mincov 80` and is dropped. This is deliberate abricate
  parity, not a rounding bug.
- **Identity is never re-filtered.** `%IDENTITY` is the BLAST `pident`, printed as `%.2f`.
  `--minid` is enforced inside blastn via `-perc_identity`.
- **Dedup, no merging.** BLAST rows with an identical `(contig, start, end)` query span are
  collapsed and the first row wins. Overlapping hits at *different* spans are all reported;
  gapit does not merge overlapping intervals. This matches abricate exactly and is a feature.
  (The opt-in `--merge-fragments` mode above is the one sanctioned exception; the default
  path never merges.)
- **Protein databases.** Databases with `dbtype prot` (for example `bacmet2`) screen through
  `blastx`, which accepts no `-perc_identity`. gapit then prints
  `--minid is not applied to protein databases (abricate parity)` on stderr and keeps going.
- **Ordering.** Within a file, rows sort by SEQUENCE then START, stable. Rows for one file are
  printed once that file finishes; files emit in input order even with `--jobs > 1`.
- **Exit codes.** 2 usage, 3 missing dependency, 4 database error, 5 input error, 1 unexpected.
  Failures print a `gapit.error/1` JSON envelope on stderr; see [./outputs.md](./outputs.md).
