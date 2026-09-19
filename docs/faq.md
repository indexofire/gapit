# FAQ

Terse answers to questions that come up in practice. Every claim is backed by
the source or a live run.

## Is gapit a drop-in abricate replacement?

For contig screening, yes: same BLAST pipeline, same hit rules, abricate-format
TSV on stdout, byte-parity checked against real abricate on a corpus. It reads
abricate-format databases (legacy `~~~` headers) as-is. The reverse does not
hold: gapit-native databases (`gapit/v1` tagged headers) cannot be read by
abricate. Details: `./databases.md`.

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

## What is the difference between `--threads` and `--jobs`?

`--threads` is BLAST worker threads inside one screening run (passed to
`-num_threads`, default 1). `--jobs` is how many input files are screened in
parallel (default 1). For N single-contig files, `--jobs N` scales; for one
huge file, `--threads` does.

## Why does reads mode reject `--format tsv`?

Reads results have no per-hit TSV semantics (they are per-gene breadth/depth
rows), so TSV/CSV is a usage error. Use `json` (default) or `md`:

```console
$ gapit screen --r1 reads.fq --read-type sr --db tinyamr --format tsv
{"schema":"gapit.error/1","code":"USAGE_ERROR","message":"--format tsv|csv is not available in reads mode (use json or md)","context":{}}
$ echo $?
2
```

Reads mode: `./reads.md`.

## What is this "MISSING_DEPENDENCY" error?

An external binary (blastn, makeblastdb, blastdbcmd, any2fasta, minimap2) is
not on PATH. gapit exits 3 with the envelope naming the binary:

```text
{"schema":"gapit.error/1","code":"MISSING_DEPENDENCY","message":"required binary not found on PATH: blastn","context":{"binary":"blastn"}}
```

Install BLAST+ and friends (pixi does this for you) and make sure they are on
PATH.

## What licenses apply to the bundled databases?

gapit itself is GPL-2.0-compatible. Database content keeps its original
upstream licenses (NCBI, CARD, CGE, and so on); gapit does not relicense it.
Details in SPEC.md §9.

## How do I add a new database to the default set?

Write a provider module, drop a snapshot archive into
`src/gapit/data/snapshots/`, and set `snapshot="<name>.tar.gz"` on the
provider so `gapit db fetch` installs offline. The header format, transforms,
and build pipeline are specified in SPEC.md §11.
