# Screening contigs

`gapit screen` aligns contig files against a gene database with BLAST and reports which genes
are present, using the same pipeline and hit rules as abricate. Output is abricate-compatible
TSV by default; JSON and Markdown are first-class alternatives (see [./outputs.md](./outputs.md)).

Database setup is covered in [./databases.md](./databases.md); a full walk-through lives in
[./quickstart.md](./quickstart.md).

## Options

Transcribed from `gapit screen --help` (gapit 0.1.0). Flags marked *reads mode* apply only when
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
| `--fofn` | path | none | File of filenames; replaces the positional FILEs. |
| `--quiet` | flag | off | Silence stderr diagnostics. |
| `--csv` | flag | off | Compat alias for `--format csv`. |
| `--noheader` | flag | off | Suppress the `#FILE ...` header row. |
| `--nopath` | flag | off | Basename the FILE column. |
| `--debug` | flag | off | Verbose stderr diagnostics; echoes each external command line. |
| `--format` | tsv\|csv\|json\|md | `tsv` | Output format (reads mode defaults to json). |
| `--r1` | str | none | *Reads mode.* Comma-separated FASTQ R1 file(s), one per lane. |
| `--r2` | str | none | *Reads mode.* Comma-separated mate FASTQ file(s); must match `--r1` count. |
| `--read-type` | sr\|map-ont\|map-hifi | `sr` | *Reads mode.* minimap2 preset. |
| `--min-breadth` | float | `90.0` | *Reads mode.* Minimum %breadth for presence. |

\* Positional FILEs or `--fofn`, or reads mode via `--r1`. Positional files and `--r1`/`--r2`
are mutually exclusive.

## Input files

any2fasta normalizes each input, so plain FASTA, gzipped and bzip2-compressed FASTA, GBK, and
EMBL all work. Gapit never loads the whole file; normalization streams into blastn. If
normalization fails (not a sequence file), gapit prints a `gapit.error/1` envelope on stderr
and exits 5.

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
    "version": "0.1.0"
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

Echoes every external command line to stderr as a `gapit: run:` line:

```console
$ gapit screen tests/data/contigs/full.fa --db tinyamr --debug 2>&1 >/dev/null
Processing: tests/data/contigs/full.fa
gapit: run: any2fasta -q -u tests/data/contigs/full.fa
gapit: run: blastn -task blastn -dust no -perc_identity 80.0 -db /tmp/gapit-demo/datadir/tinyamr/sequences -outfmt '6 qseqid qstart qend qlen sseqid sstart send slen sstrand evalue length pident gaps gapopen stitle' -num_threads 1 -evalue 1E-20 -culling_limit 1 -max_target_seqs 10000
Found 1 genes in tests/data/contigs/full.fa
```

stdout is unchanged by `--debug`; it is a stderr-only flag.

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
- **Protein databases.** Databases with `dbtype prot` (for example `bacmet2`) screen through
  `blastx`, which accepts no `-perc_identity`. gapit then prints
  `--minid is not applied to protein databases (abricate parity)` on stderr and keeps going.
- **Ordering.** Within a file, rows sort by SEQUENCE then START, stable. Rows for one file are
  printed once that file finishes; files emit in input order even with `--jobs > 1`.
- **Exit codes.** 2 usage, 3 missing dependency, 4 database error, 5 input error, 1 unexpected.
  Failures print a `gapit.error/1` JSON envelope on stderr; see [./outputs.md](./outputs.md).
