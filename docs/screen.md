# Screening contigs

`gapit screen` aligns contig files against a gene database with BLAST and reports which genes
are present, using the same pipeline and hit rules as abricate. Output is abricate-compatible
TSV by default; JSON and Markdown are first-class alternatives (see [./outputs.md](./outputs.md)).

Database setup is covered in [./databases.md](./databases.md); a full walk-through lives in
[./quickstart.md](./quickstart.md).

## Options

Transcribed from `gapit screen --help` (gapit 0.5.2). Flags marked *reads mode* apply only when
you pass `--r1`/`--r2`; they are documented in [./reads.md](./reads.md).

Every option also accepts the single-dash short form listed in the **Short** column (e.g.
`-d` for `--db`); long forms remain the canonical spelling, and boolean negative halves like
`--no-merge-fragments` stay long-only. `--db` alone carries a second, word-style short alias:
`-db ncbi` (space-separated; an attached `-dbncbi` still binds to `-d` plus a value). A
positional argument that names a known database but is missing on disk is rejected with a
targeted usage error pointing at `--db`.

| Flag | Short | Type | Default | Description |
|---|---|---|---|---|
| `FILE...` | — | path(s) | required* | Input contig file(s) to screen. An all-FASTQ wildcard enters reads mode with samples auto-paired from filenames (see [./reads.md](./reads.md)). |
| `--db` | `-d`, `-db` | str | required | Database to screen against (datadir subdir). No default: a screen never silently materializes a bundled database — pick one explicitly (`gapit db list`). |
| `--datadir` | `-D` | path | `$GAPIT_DATADIR`, then `~/.local/share/gapit/db` | Database directory. |
| `--minid` | `-i` | float | `80.0` | Minimum %identity, `0 < x <= 100`. Enforced inside BLAST via `-perc_identity`. |
| `--mincov` | `-c` | float | `80.0` | Minimum %coverage, `0 <= x <= 100`. Post-filter on the unrounded float. |
| `--threads` | `-t` | int | `1` | BLAST worker threads (passed to `-num_threads`). |
| `--jobs` | `-j` | int | `1` | Screen N input files concurrently. Output order is always input order. |
| `--merge-fragments` | `-m` | flag | off | Merge gene fragments split across contigs (gapit extension; see below). |
| `--fofn` | `-F` | path | none | File of filenames; replaces the positional FILEs. |
| `--quiet` | `-q` | flag | off | Silence stderr diagnostics. |
| `--noheader` | `-n` | flag | off | Suppress the `#FILE ...` header row. |
| `--nopath` | `-p` | flag | off | Basename the FILE column. |
| `--debug` | `-v` | flag | off | Verbose stderr diagnostics; echoes each external command line. |
| `--format` | `-f` | tsv\|csv\|json\|md | `tsv` | Output format (tsv is the default everywhere: reads mode streams the table per completed file or sample). json/md are the explicit agent opt-ins; json is written once at the end (single document). |
| `--output` | `-o` | path | stdout | Write the report to PATH instead of stdout (truncates any existing file). Streaming formats flush per file; stdout then carries no data. |
| `--aligner` | `-a` | blastn\|minimap2 | input-based | Alignment engine (default: blastn for contig files, minimap2 for `--r1`/`--r2` reads). `--aligner minimap2` routes positional FASTA assemblies through the minimap2 engine (FASTA content required; see [./reads.md](./reads.md)). |
| `--r1` | `-1` | str | none | *Reads mode.* Reads or assembly FASTA file(s), comma-separated, one per lane. |
| `--r2` | `-2` | str | none | *Reads mode.* Comma-separated mate FASTQ file(s); must match `--r1` count. |
| `--read-type` | `-x` | sr\|map-ont\|map-hifi | `sr` for FASTQ, `map-ont` for FASTA | *Reads mode.* minimap2 preset; resolved from the detected input when omitted. |
| `--min-breadth` | `-b` | float | `90.0` | *Reads mode.* Minimum %breadth for presence. |
| `--min-gene-cov` | `-g` | float | `90.0` | *Cluster databases only.* Minimum %coverage for a gene `present` verdict. |
| `--min-gene-id` | `-G` | float | `90.0` | *Cluster databases only.* Minimum %identity for a gene `present` verdict. |
| `--min-cluster-cov` | `-C` | float | `96.0` | *Cluster databases only.* Minimum locus %coverage for a best-locus call. |

\* Positional FILEs or `--fofn`, or reads mode via `--r1`. Positional files and `--r1`/`--r2`
are mutually exclusive.

\* Positional FILEs or `--fofn`, or reads mode via `--r1`. Positional files and `--r1`/`--r2`
are mutually exclusive.

## Input files

Normalization is native (no external `any2fasta`): plain FASTA, gzipped and bzip2-compressed
FASTA, FASTQ, GBK, and EMBL all work. The converted FASTA is buffered in memory (genome-scale
assemblies are a few MB) and fed to blastn on stdin. If normalization fails (not a sequence
file), gapit prints a `gapit.error/1` envelope on stderr and exits 5.

One exception predates normalization: when **every** positional file is FASTQ (`.fastq`/`.fq`
± `.gz`, or content-sniffed FASTQ under an ambiguous extension) and no `--aligner` is given,
the whole invocation routes to the reads engine — `gapit screen -d ecoli_dec *.gz` screens
the glob as auto-paired samples instead of contigs. Mixing FASTA and FASTQ positionals is a
usage error naming the reads files. See [./reads.md](./reads.md) for the pairing
conventions and the per-sample output.

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
    "version": "0.5.2"
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

### Typed gene databases: designation is a second command

When `--db` names a gene database built with a `typing.json` (the bundled `ecoli_dec`, or
any `gapit db build --typing` build), screening is pure gene detection: the output is the
frozen 15-column abricate table byte-identically in every format — typed and untyped gene
databases are indistinguishable on the screen surface. Designations come from the
two-stage pipeline's second stage: write the table with `--output` (or redirect), then run
[`gapit typing`](./typing.md#the-two-stage-designation-workflow) on it — or pipe it
straight in, the canonical one-liner:

```console
$ gapit screen dec_s3_stx2a_escV_aggR_uidA.fasta dec_s2_pic_astA_uidA.fasta --db ecoli_dec --output dec.tsv --nopath --quiet
$ gapit typing dec.tsv --quiet
FILE	SCHEME	PHENOTYPE	CONFIDENCE	SCORE	RUNNER_UP	NOTES	GENES
dec_s3_stx2a_escV_aggR_uidA.fasta	gb4789_6	EHEC	high	1.0000	EAEC (1.0000)	GB 4789.6-2016: EHEC = stx (any subunit) + escV (LEE); severity order EHEC>STEC/EPEC>ETEC>EIEC>EAEC: rules are declared in severity order so hybrids surface as runner_up	aggR;escV;stx2a;uidA
dec_s3_stx2a_escV_aggR_uidA.fasta	risk_monitoring	EHEC	high	1.0000	EAEC (1.0000)	GB 4789.6-2016: EHEC = stx (any subunit) + escV (LEE); severity order EHEC>STEC/EPEC>ETEC>EIEC>EAEC: rules are declared in severity order so hybrids surface as runner_up	aggR;escV;stx2a;uidA
dec_s2_pic_astA_uidA.fasta	gb4789_6	EAEC	high	1.0000	EHEC (0.0000)	GB 4789.6: any of aggR/pic/astA	astA;pic;uidA
dec_s2_pic_astA_uidA.fasta	risk_monitoring	non-DEC	low	0.0000	STEC (0.0000)	GB 4789.6-2016: EHEC = stx (any subunit) + escV (LEE); severity order EHEC>STEC/EPEC>ETEC>EIEC>EAEC: rules are declared in severity order so hybrids surface as runner_up	astA;pic;uidA
```

The piped form skips the intermediate file entirely — screen's stdout is typing's stdin
(stdout purity holds on both sides: data on the pipe, diagnostics on stderr):

```console
$ gapit screen dec_s3_stx2a_escV_aggR_uidA.fasta --db ecoli_dec --nopath --quiet | gapit typing --quiet
```

`gapit typing` reads stdin whenever it is not a terminal; `gapit typing -` is the explicit
stdin marker (valid with a terminal attached — it reads until EOF). The full input
contract lives in the [typing page](./typing.md#the-two-stage-designation-workflow).

The same pipe feeds the matrix view. Screening many assemblies in one call batched all
their rows into one table (every row carries the same `DATABASE`; the `FILE` column spans
the whole batch), so [`gapit summary`](./summary.md) summarizes it straight off the pipe
into the gene×file presence matrix:

```console
$ gapit screen -d ecoli_dec *.fna --quiet | gapit summary
```

Cluster databases are the exception: their typing is integrated into the screen itself
(the typed cluster TSV's `PHENOTYPE` column, [below](#cluster-databases-kind-cluster)).

## Streaming and `--output`

Long batch runs give feedback as they go, and the report can go to a file:

- **tsv/csv/md stream per file.** The header (or Markdown frontmatter) prints once screening
  starts, and each file's rows/section print the moment that file finishes — in input order,
  always. Under `--jobs N > 1` emission is head-of-line: file *i*'s output waits until files
  1..*i* are all done (the pool yields in input order), so bytes on stdout are identical to
  the sequential run. The positional-FASTQ wildcard path (`screen -d db *.fastq.gz`,
  [./reads.md](./reads.md)) has the same contract sample-wise: `--jobs` parallelizes
  samples, md streams the static frontmatter plus one `## <sample>` section per completed
  sample, and `--output` persists the streamed prefix on a mid-batch failure.
- **md streams; json is written once, at the end.** The Markdown frontmatter is STATIC
  metadata (schema, tool, `created_at`, db, thresholds — no run totals), so it can lead the
  document and every file's section follows the moment the file completes, exactly like the
  tsv rows. On cluster databases each file's section carries its own summary row and gene
  table. Run totals (`files`, `hits`, ...) live in the JSON document, which is a single
  object written once at the end.
- **`--output PATH`** (all engines, including reads mode) writes the report to PATH instead of
  stdout: the file opens on the first output byte (truncating any existing file — v1 overwrite
  semantics, never append), every streamed chunk is flushed, and stdout then carries **no
  data** (stderr diagnostics are unchanged).
- **Errors mid-batch.** If file *k* fails after files 1..*k-1* streamed, the already-emitted
  output persists — on stdout it is already printed; with `--output` the file keeps the
  header/frontmatter plus files 1..*k-1* — and then the typed `gapit.error/1` envelope prints
  on stderr with the documented exit code. A run that fails before any output (usage,
  dependency, db errors) creates no `--output` file at all.

### `--output`

```console
$ gapit screen tests/data/contigs/full.fa tests/data/contigs/gap.fa --db tinyamr --output report.tsv
Processing: tests/data/contigs/full.fa
Found 1 genes in tests/data/contigs/full.fa
Processing: tests/data/contigs/gap.fa
Found 1 genes in tests/data/contigs/gap.fa
$ cat report.tsv
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
tests/data/contigs/full.fa	contig1	1	79	+	tetA	1-79/79	===============	0/0	100.00	100.00	tinyamr	NC_000913.3:100-900	tetracycline efflux pump TetA	TETRACYCLINE
tests/data/contigs/gap.fa	contig1	1	97	+	sul1	1-94/94	========/======	1/3	100.00	96.91	tinyamr	U12338.4:1-940	sulfonamide-resistant dihydropteroate synthase Sul1	SULFONAMIDE
```

The `Processing:`/`Found` lines are stderr; stdout is empty and `report.tsv` holds exactly the
bytes the same run would print.

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
  printed (or flushed to `--output`) once that file finishes; files emit in input order even
  with `--jobs > 1` (head-of-line: file i waits for files 1..i).
- **Exit codes.** 2 usage, 3 missing dependency, 4 database error, 5 input error, 1 unexpected.
  Failures print a `gapit.error/1` JSON envelope on stderr; see [./outputs.md](./outputs.md).
