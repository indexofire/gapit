# FAQ

Terse answers to questions that come up in practice. Every claim is backed by
the source or a live run.

## Is gapit a drop-in abricate replacement?

Contig screening keeps the abricate surface on purpose: same BLAST pipeline, same hit
rules, abricate-format TSV on stdout, byte-checked against real abricate on a corpus.
From v0.5.0 that baseline is frozen rather than a release gate — gapit pursues its own
contract (native JSON/Markdown, reads and cluster engines, typing), and the parity
harness (`pixi run -e parity parity`) serves as a regression reference for the frozen
surface. gapit reads abricate-format databases (legacy `~~~` headers) as-is. The
reverse does not hold: gapit-native databases (`gapit/v1` tagged headers) cannot be
read by abricate. Details: `./databases.md`.

## Why does a hit showing 80.00% coverage get filtered?

The `--mincov 80` comparison runs on the unrounded float; the TSV column is
display-rounded. A real case: a 1000 nt gene with a 1 nt gap in a 25000 nt
reference gives 100*(20000-1)/25000 = 79.996%, which prints as `80.00` but is
dropped. This matches abricate exactly. Background: `./screen.md`.

## Why do I get zero hits?

Check the thresholds. Defaults are `--minid 80` (enforced inside blastn via
`-perc_identity`) and `--mincov 80` (post-filter on the unrounded coverage,
see the question above). Lower them for divergent or partial genes. Note that
an input with no matches is not an error: the header row still prints and the
exit code is 0.

```console
$ gapit screen none.fa --db tinyamr
Processing: none.fa
Found 0 genes in none.fa
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
$ echo $?
0
```

## How do I update a database?

Refetch over the old one with `--force`:

```console
$ gapit db fetch ncbi --force
```

Without `--force`, fetching into an existing database fails instead of
overwriting.

## Where is my datadir?

Resolution order: `--datadir` on the command line, then `$GAPIT_DATADIR`, then
`~/.local/share/gapit/db`. A resolved path that does not exist is an error
(`DATADIR_NOT_FOUND`, exit 4). See `./databases.md`.

## Can I use my existing abricate databases?

Yes, unchanged: point gapit at the datadir that abricate built and screen.
If the BLAST indices are missing or stale, rebuild them all with:

```console
$ gapit setupdb
```

Only gapit-native output goes the other way (abricate cannot read `gapit/v1`
databases).

## How do I screen gene clusters or Kaptive loci?

Cluster databases are a second database kind, built from GenBank/GFF files
with `gapit db build NAME loci.gbk` (optionally `--typing FILE` for phenotype
calls) or fetched from the seven kaptive providers (`kpsc_k`, `kpsc_o`, `kosc_k`,
`kosc_o`, `ab_k`, `ab_o`, `ecoli_kps` — the official Kaptive v3 install keywords).
`gapit screen assembly.fa --db kpsc_k` dispatches to the
minimap2 cluster engine and reports one best-locus call per file, with
per-gene verdicts and, on typed databases, a phenotype. Details:
`./screen.md` and `./databases.md`.

## How do I designate a pathotype or serogroup?

Designation is two-stage on gene databases: screen with a typed database, then read
the table back. The bundled `ecoli_dec` panel carries the dual-scheme DEC designation:
`gapit screen -d ecoli_dec sample.fna -o result.tsv`, then `gapit typing result.tsv`
— or pipe one straight into the other. Cluster databases are the exception: their
designation is integrated into the screen itself (the `PHENOTYPE` column). Details:
`./typing.md`.

## What is the difference between `--threads` and `--jobs`?

`--threads` is BLAST worker threads inside one screening run (passed to
`-num_threads`, default 1). `--jobs` is how many input files are screened in
parallel (default 1). For N single-contig files, `--jobs N` scales; for one
huge file, `--threads` does.

## What format does reads mode output?

TSV by default — tsv is the human default on every gapit surface, and json/md are the
agent opt-ins. Reads results are per-gene breadth/depth rows, so the streaming table
puts the sample key in the first `#SAMPLE` column; `--format json` (or `md`) gives the
versioned `gapit.reads/1` document instead:

```console
$ gapit screen --r1 tests/data/reads/tetx_full.fq --read-type sr --db tinyreads
#SAMPLE	GENE	BREADTH%	DEPTH	READS	PRESENT	DATABASE	ACCESSION	PRODUCT
tests/data/reads/tetx_full.fq	tetX	100.00	2.30	12	yes	tinyreads	SYN-001	extended resistance determinant tetX
$ gapit screen --r1 tests/data/reads/tetx_full.fq --read-type sr --db tinyreads --format json
{"schema": "gapit.reads/1", ...}
```

Reads mode: `./reads.md`.

## Can I screen a wildcard of FASTQ files?

Yes: when every positional file is FASTQ, `gapit screen` enters reads mode and
auto-pairs samples from the filenames (`_R1`/`_R2`, `_1`/`_2` conventions), so
`gapit screen -d ecoli_dec *.fq.gz -j 4` screens all samples concurrently with one
command. Mixing FASTA and FASTQ positionals is a usage error; the pairing conventions
and the per-sample output are documented in `./reads.md`.

## What are per-gene identity floors?

A database-side `floors.json` (`gapit.floors/1`) declares a minimum alignment identity
per gene for reads-mode presence: alignments below the floor are dropped before
breadth and depth are aggregated, so a close-but-not-exact homolog no longer over-calls
the gene. The real case that motivated it: the SPATE-homolog `pic` in the bundled
`ecoli_dec` panel read about 86.6% identity at 97.6% breadth and crossed the 90%
presence threshold, while the blastn contig path correctly rejected it; the shipped
floor (`pic` at 90) makes reads mode agree. Install one with `gapit db build --floors
FILE` or by dropping the file into the db directory; a database without one screens
byte-identically to before. Details: `./reads.md` and `./databases.md`.

## What is this "MISSING_DEPENDENCY" error?

An external binary (blastn, makeblastdb, blastdbcmd, minimap2) is
not on PATH. gapit exits 3 with the envelope naming the binary:

```text
{"schema":"gapit.error/1","code":"MISSING_DEPENDENCY","message":"required binary not found on PATH: blastn","context":{"binary":"blastn"}}
```

Install BLAST+ and friends (pixi does this for you) and make sure they are on
PATH.

## What licenses apply to the databases?

gapit itself is MIT-licensed. Six audited, permissively licensed bundles ship inside
the wheel (the public-domain `ecoli_dec` and `lm_doumith` panels plus the `ncbi`,
`resfinder`, `ecoh`, `upec_expec_vf` snapshots: Apache-2.0 / BSD-3-Clause / MIT) and
materialize into the datadir on first use; every other provider downloads from upstream
at fetch time, and the content keeps its original
licenses (CARD's McMaster non-commercial terms, VFDB's CC BY-NC,
CGE, Kaptive's GPL-3.0, and so on); gapit does not relicense it. Providers that pin a
license expose it in `gapit db list --json`. Details in SPEC.md §9 and the
[bundled section](./databases.md#bundled-databases-install-time-ready).

## How do I add a new database to the default set?

Write a provider module and fetch it by name: every provider downloads from its
upstream source at fetch time, transforms the records, and builds locally. Nothing
beyond the six audited bundles ships in the wheel (upstream licenses such as CARD's or
VFDB's forbid redistribution inside an MIT-licensed distribution). The header format,
transforms, and build pipeline are specified in SPEC.md §11.
