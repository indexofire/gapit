# Databases

A gapit database is a directory of reference genes you screen contigs and reads against.
This page covers where databases live, how to get them, what gapit writes on disk, and how
to bring your own. For how screening uses them, see [screen.md](./screen.md).

## Where databases live

Every database sits in a **datadir**, one subdirectory per database. gapit resolves the
datadir in this order:

1. `--datadir` on the individual call
2. `$GAPIT_DATADIR`
3. `~/.local/share/gapit/db`

The resolved path must exist. If it doesn't, commands fail with a typed error and exit 4:

```console
$ gapit db list --datadir /no/such/dir
{"schema":"gapit.error/1","code":"DATADIR_NOT_FOUND","message":"datadir does not exist: /no/such/dir","context":{"datadir":"/no/such/dir"}}
```

## Listing what's installed

`gapit db list` shows every known provider and whether it is installed:

```console
$ gapit db list
PROVIDER	STATUS	DBTYPE	DESCRIPTION
argannot	installed (2224)	nucl	ARG-ANNOT acquired resistance genes
bacmet2	installed (746)	prot	BacMet2 experimentally confirmed biocide/resistance genes (protein)
card	installed (6059)	nucl	CARD protein homolog resistance models
ecoh	installed (597)	nucl	E. coli O and H antigens (srst2 EcOH)
ecoli_vf	installed (2701)	nucl	E. coli virulence factors (phac-nml)
megares	installed (7425)	nucl	MEGARes antimicrobial resistance genes
ncbi	installed (8373)	nucl	NCBI AMRFinderPlus (reference finder) curated AMR
plasmidfinder	installed (488)	nucl	CGE PlasmidFinder replicons
resfinder	installed (3206)	nucl	CGE ResFinder acquired resistance genes
upec_expec_vf	installed (77)	nucl	UPEC/ExPEC virulence genes (FordeGenomics)
vfdb	installed (4769)	nucl	VFDB virulence factors (set A, nucleotide)
victors	installed (4402)	nucl	Victors virulence factors
```

STATUS reads `installed (N)` when `<datadir>/<name>/gapit-manifest.json` exists, with N the
record count, otherwise `available`. For agents, `--json` emits a `gapit.dblist/1` document
(first three of twelve providers shown, output trimmed):

```console
$ gapit db list --json
{
  "schema": "gapit.dblist/1",
  "providers": [
    {
      "name": "argannot",
      "description": "ARG-ANNOT acquired resistance genes",
      "dbtype": "nucl",
      "installed": true,
      "records": 2224
    },
    {
      "name": "bacmet2",
      "description": "BacMet2 experimentally confirmed biocide/resistance genes (protein)",
      "dbtype": "prot",
      "installed": true,
      "records": 746
    },
    {
      "name": "card",
      "description": "CARD protein homolog resistance models",
      "dbtype": "nucl",
      "installed": true,
      "records": 6059
    }
  ]
}
```

| Field | Type | Meaning |
|---|---|---|
| `schema` | string | Always `gapit.dblist/1` |
| `providers` | array | One entry per provider |
| `providers[].name` | string | Provider name, as passed to `gapit db fetch` |
| `providers[].description` | string | Short content summary |
| `providers[].dbtype` | string | `nucl` (screened with blastn) or `prot` (screened with blastx) |
| `providers[].installed` | boolean | True when a manifest exists in the datadir |
| `providers[].records` | integer | Record count, omitted when the database isn't installed |

`gapit list` gives the abricate-compatible view of installed databases; `gapit list --json`
returns a `gapit.list/1` document. See [outputs.md](./outputs.md).

## Checking database freshness

`gapit db outdated` reports every installed database's age and flags two update conditions.
Against a fully installed datadir (output trimmed to three of twelve rows):

```console
$ gapit db outdated --days 30
NAME	FETCHED_AT	AGE_DAYS	STATUS
argannot	2026-09-17T23:22:06Z	2.58	ok
bacmet2	2026-09-17T23:13:51Z	2.58	ok
card	2026-09-17T23:14:25Z	2.58	ok
...
```

A database installed long ago and superseded by the bundled snapshot reports both flags:

```console
$ gapit db outdated --datadir /tmp/opencode/gapit-outdated-demo
NAME	FETCHED_AT	AGE_DAYS	STATUS
card	2020-01-01T00:00:00Z	2454.55	stale+snapshot-update
```

| Status | Meaning |
|---|---|
| `ok` | Fresh enough and no newer bundle |
| `stale` | `age_days` past `--days` (default 90; `--days 0` marks everything stale) |
| `snapshot-update` | The provider's bundled snapshot is newer than the installed copy (card, vfdb) |
| `stale+snapshot-update` | Both of the above |

Staleness is a report, never an error state: the command exits 0 however stale things are.
Exit 4 (`DATADIR_NOT_FOUND`, `DATADIR_EMPTY`) covers a missing datadir or one with no installed
databases; an unparseable `fetched_at` is `MANIFEST_MALFORMED` (exit 5). For agents, `--json`
emits a `gapit.dboutdated/1` document (a CLI listing like `gapit.dblist/1`, not registered with
`gapit schema`):

```console
$ gapit db outdated --days 3 --json
{
  "schema": "gapit.dboutdated/1",
  "databases": [
    {
      "db": "argannot",
      "fetched_at": "2026-09-17T23:22:06Z",
      "age_days": 2.58,
      "status": "ok",
      "upstream_version": ""
    },
    {
      "db": "card",
      "fetched_at": "2026-09-17T23:14:25Z",
      "age_days": 2.58,
      "status": "ok",
      "upstream_version": ""
    }
  ]
}
```

## Searching records across databases

`gapit db search TERM` looks up genes in the `records.jsonl` truth store of every installed
database — no BLAST, just a fast scan. Matching is case-insensitive substring by default;
`--exact` switches to full-field equality:

```console
$ gapit db search ctx-m --limit 3
argannot	(Bla)blaCTX-M-1	X92506:63-938		(Bla)blaCTX-M-1	876
argannot	(Bla)blaCTX-M-10	AF255298:1-873		(Bla)blaCTX-M-10	873
argannot	(Bla)blaCTX-M-100	FR682582:1-876		(Bla)blaCTX-M-100	876
$ gapit db search "blaCTX-M-1" --field gene --exact
ncbi	blaCTX-M-1	NG_048897.1	CEPHALOSPORIN	extended-spectrum class A beta-lactamase CTX-M-1	876
$ gapit db search virulence --db vfdb --field function --limit 2
vfdb	AAA92657	AAA92657	virulence	(AAA92657) unknown protein [TraJ (VF0241) - Invasion (VFC0083)] [Escherichia coli]	606
vfdb	AAC38364	AAC38364	virulence	(AAC38364) Orf1 [Ler (VF0189) - Regulation (VFC0301)] [Escherichia coli O127:H6 str. E2348/69]	390
```

Rows are `DB\tGENE\tACCESSION\tFUNCTION\tPRODUCT\tLENGTH`, streamed in database-then-file
order; `FUNCTION` joins the record's function classes with `;`. `--field` picks
`gene|accession|function|product|any` (default `any` searches all of them, one function class
at a time). `--limit N` caps the output (default 100; `0` = unlimited) and stderr notes a
truncation (`--quiet` silences it); zero hits exit 0 with empty stdout. `--json` prints one
JSON object per hit with the same fields as snake_case (`function` as an array):

```console
$ gapit db search "tet(M)" --field gene --exact --json | head -1
{"db":"card","gene":"tet(M)","accession":"AB039845.1:25-1945","function":["tetracycline"],"product":"Tet(M) is a ribosomal protection protein that confers tetracycline resistance. It is found on transposable DNA elements and its horizontal transfer between bacterial species has been documented.","length":1920}
```

An installed database whose `records.jsonl` is missing is skipped with a stderr warning during
a full scan, but targeting it explicitly (`--db NAME`) fails with exit 4 `DB_INCOMPLETE`; an
unknown `--db NAME` is a usage error (exit 2) listing what is installed.

## Providers

Twelve providers ship with gapit. `card` and `vfdb` also ship as bundled snapshots inside
the package, so they install with zero network access.

| Name | Content | dbtype |
|---|---|---|
| `ncbi` | NCBI AMRFinderPlus (reference finder) curated AMR (default db) | nucl |
| `card` | CARD protein homolog resistance models | nucl |
| `resfinder` | CGE ResFinder acquired resistance genes | nucl |
| `argannot` | ARG-ANNOT acquired resistance genes | nucl |
| `plasmidfinder` | CGE PlasmidFinder replicons | nucl |
| `megares` | MEGARes antimicrobial resistance genes | nucl |
| `ecoh` | E. coli O and H antigens (srst2 EcOH) | nucl |
| `vfdb` | VFDB virulence factors (set A, nucleotide) | nucl |
| `ecoli_vf` | E. coli virulence factors (phac-nml) | nucl |
| `bacmet2` | BacMet2 experimentally confirmed biocide/resistance genes (protein) | prot |
| `victors` | Victors virulence factors | nucl |
| `upec_expec_vf` | UPEC/ExPEC virulence genes (FordeGenomics) | nucl |

Protein databases (`bacmet2`) screen through `blastx`; nucleotide ones through `blastn`.

## Fetching databases

### Bundled snapshots, zero network

Bare `gapit db fetch` (no name) installs the default set, `card` then `vfdb`, from the
snapshots bundled in the package. Nothing touches the network: the snapshot archive carries
`records.jsonl` plus the manifest, and gapit rebuilds `sequences` and the BLAST index
locally. That rebuild is deterministic and fast, and it matches the BLAST
version actually installed on your machine.

The bundled data does not rot: a scheduled workflow
(`.github/workflows/snapshot-refresh.yml`) re-fetches card and vfdb from upstream monthly
and opens a pull request whenever the records changed. That PR is the review gate — a human
signs off on the data update before the new tars merge. Outside GitHub Actions,
`gapit db fetch <name> --from-source` remains the manual upstream path.

### Fetching a named provider

`gapit db fetch <name>` runs the full pipeline for one provider. Here is `card` from its
bundled snapshot (stderr progress lines, then a one-line JSON receipt on stdout):

```console
$ gapit db fetch --datadir /tmp/opencode/gapit-dbs-demo card
gapit: installed card from bundled snapshot card.tar.gz
gapit: generated /tmp/opencode/gapit-dbs-demo/card/sequences
gapit: self-check passed for card
gapit: BLAST index built (nucl)
{"db":"card","records":6059,"dbtype":"nucl","destination":"/tmp/opencode/gapit-dbs-demo/card"}
```

| Receipt field | Type | Meaning |
|---|---|---|
| `db` | string | Provider name |
| `records` | integer | Records written to `records.jsonl` |
| `dbtype` | string | `nucl` or `prot` |
| `destination` | string | Installed database directory |

Providers without a bundled snapshot download from their upstream source at fetch time, so
they need network access. Unknown names are rejected before the datadir is even resolved:

```console
$ gapit db fetch nosuchdb
{"schema":"gapit.error/1","code":"USAGE_ERROR","message":"unknown provider: nosuchdb (available: argannot, bacmet2, card, ecoh, ecoli_vf, megares, ncbi, plasmidfinder, resfinder, upec_expec_vf, vfdb, victors)","context":{"provider":"nosuchdb"}}
```

### Refetching and forcing upstream

An existing database directory is never overwritten silently. Refetching without flags
fails with exit 4:

```console
$ gapit db fetch --datadir /tmp/opencode/gapit-dbs-demo card
{"schema":"gapit.error/1","code":"DB_ALREADY_EXISTS","message":"won't overwrite existing database card (use --force)","context":{"db":"card"}}
```

`--force` deletes and rebuilds in place:

```console
$ gapit db fetch --datadir /tmp/opencode/gapit-dbs-demo --force card
{"db":"card","records":6059,"dbtype":"nucl","destination":"/tmp/opencode/gapit-dbs-demo/card"}
```

`--from-source` skips the bundled snapshot and forces the upstream download even when a
snapshot exists. Conversely, if a provider has no snapshot archive, the fetch falls back to
the network path on its own.

## Installing a local file

`gapit db install SOURCE --sha256 HASH --output TARGET` is a verified local-file install.
It knows nothing about providers, archives, or sequence content: it streams SOURCE through
SHA256, compares the digest against `--sha256`, and only then atomically replaces TARGET.
A failed check leaves any existing TARGET untouched.

```console
$ sha256sum src/gapit/data/snapshots/card.tar.gz
65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025  src/gapit/data/snapshots/card.tar.gz
$ gapit db install src/gapit/data/snapshots/card.tar.gz \
    --sha256 65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025 \
    --output /tmp/opencode/card-copy.tar.gz
{"destination":"/tmp/opencode/card-copy.tar.gz","sha256":"65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025"}
```

A wrong digest aborts with exit 5:

```console
$ gapit db install src/gapit/data/snapshots/card.tar.gz \
    --sha256 0000000000000000000000000000000000000000000000000000000000000000 \
    --output /tmp/opencode/card-bad.tar.gz
{"schema":"gapit.error/1","code":"CHECKSUM_MISMATCH","message":"SHA256 mismatch for src/gapit/data/snapshots/card.tar.gz","context":{"source":"src/gapit/data/snapshots/card.tar.gz","expected":"0000000000000000000000000000000000000000000000000000000000000000","actual":"65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025"}}
```

## What's inside a gapit-built database

```
<datadir>/<name>/
  records.jsonl         truth source: one Record JSON object per line
  sequences             generated FASTA projection (gapit/v1 tagged headers)
  sequences.n*|p*       BLAST index built from sequences
  gapit-manifest.json   provenance sidecar, written last
```

### records.jsonl

This is the editable truth. Every screening artifact is generated from it, so editing
records and refetching with `--force` regenerates everything downstream. One JSON object
per line; the first card record, verbatim:

```json
{"db":"card","gene":"23S_rRNA_(adenine(2058)-N(6))-methyltransferase_Erm(A)","sequence":"ATGAAACAGAAAAACCCGAAAAATACGCAAAATTTCATTACATCTAAAAAGCATGTAAAGGAAATATTAAAATATACGAATATCAATAAACAAGATAAAATAATAGAAATTGGGTCAGGAAAAGGACATTTTACCAAGGAACTTGTGGAAATGAGTCAACGGGTGAATGCTATAGAGATTGATGAAGGTTTATGTCATGCCACGAAAAAAGCAGTTGAACCTTTTCAGAATATAAAAGTTATTCATGAGGATATTTTGAAGTTTAGCTTTCCTAAAAATACAGACTATAAAATATTTGGTAATATTCCCTACAATATTAGTACTGATATTGTAAAAAAGATTGCTTTTGATAGTCAAGCGAAATATAGCTACCTTATTGTAGAGAGGGGATTTGCTAAAAGGTTGCAAAATACCCAACGAGCTTTAGGTTTGCTGTTAATGGTGGAAATGGATATAAAAATTCTTAAAAAAGTGCCACGAGCATATTTTCACCCTAAGCCTAATGTAGATTCTGTATTGATTGTACTTGAAAGGCATAAACCATTTATTTTAAAGAAGGACTACAAAAAGTATAGATTTTTCGTTTATAAATGGGTAAACAGGGAATATCATGTTCTTTTTACTAAAAATCAATTAAGACAGGTGCTGAAGCATGCGAATGTTACTGATCTTGATAAATTATCCAATGAACAATTTTTGTCTGTTTTCAATAGTTACAAATTATTTCAATAA","accession":"AF002716.1:210-942","function":["lincosamide","macrolide","streptogramin"],"product":"Variant of ErmA (ARO:3000347) found in Streptococcus pyogenes. Confers the MLSb phenotype.","source_id":"3005099"}
```

| Field | Type | Meaning |
|---|---|---|
| `db` | string | Source database name |
| `gene` | string | Gene name, the identity gapit reports on |
| `sequence` | string | Normalized sequence (uppercase; `N`/`X` for ambiguous bases) |
| `accession` | string | Accession or coordinate span, empty when unpublished |
| `function` | array of string | Functional categories, sorted (see the vocabulary below) |
| `product` | string | Product description, `n/a` when absent |
| `source_id` | string | Provider-native record identifier |

### gapit-manifest.json

The manifest certifies the build: what was fetched, when, from where, and with which tool
versions. A real one, from the plasmidfinder database:

```json
{
  "schema": "gapit.manifest/1",
  "name": "plasmidfinder",
  "source_urls": [
    "https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/HEAD.zip"
  ],
  "fetched_at": "2026-09-17T23:11:36Z",
  "sha256": "a18dc9dab762fe1e23a90c11314ee179c37c59dad833a5f465a78dc245970c2a",
  "n_records": 488,
  "dbtype": "nucl",
  "header_format": "gapit/v1",
  "upstream_version": "",
  "tool": {
    "name": "gapit",
    "version": "0.2.0"
  },
  "makeblastdb_version": "blastn: 2.17.0+",
  "minimap2_version": "2.31-r1302"
}
```

| Field | Type | Meaning |
|---|---|---|
| `schema` | string | Always `gapit.manifest/1` |
| `name` | string | Database name |
| `source_urls` | array of string | Upstream URLs the content came from |
| `fetched_at` | string | ISO-8601 UTC timestamp, second precision |
| `sha256` | string | Digest of the built `sequences` file, 64 lowercase hex chars |
| `n_records` | integer | Record count, matches `records.jsonl` |
| `dbtype` | string | `nucl` or `prot`, chosen explicitly at build time |
| `header_format` | string | Always `gapit/v1` for gapit-built databases |
| `upstream_version` | string | Upstream release label when the source has one |
| `tool` | object | `{name, version}` of the gapit that built the database |
| `makeblastdb_version` | string | BLAST+ version that built the index |
| `minimap2_version` | string | minimap2 version in the build environment (reads mode indexes in memory; no `.mmi` is built) |

`records.jsonl` and the manifest are file contracts. They never appear on stdout and are
not registered with `gapit schema`.

## Header formats and compatibility

gapit reads two header formats and detects them per record:

- **abricate legacy**: `>DB~~~GENE~~~ACCESSION~~~RESISTANCE PRODUCT`. Fields may be
  missing; missing trailing fields decode to empty strings.
- **gapit-native `gapit/v1`**: `>gapit|db=<v>|gene=<v>|acc=<v>|func=<v> PRODUCT` with
  percent-encoded values, so gene names and products survive BLAST and minimap2 byte for
  byte, spaces and pipes included.

Both decode into the same four identity fields, and the TSV/JSON outputs are identical in
shape either way. For native databases the `RESISTANCE` column carries the `func` value.

> **Warning:** abricate cannot read gapit-built databases. The `gapit/v1` header format is
> not the `~~~` convention. gapit reading abricate datadirs is one-directional by design;
> keep a legacy copy if a downstream tool still runs abricate.

### The `func` vocabulary

Each native provider fills `function` from a fixed vocabulary, and these values flow
through to the TSV `RESISTANCE` column and the JSON `resistance` field:

| Provider(s) | `func` value(s) |
|---|---|
| `ncbi`, `resfinder`, `argannot`, `card` | The source's antibiotic classes |
| `megares` | The MEGARes class field |
| `ecoh` | `H-antigen`, `O-antigen`, or an antigen fallback derived from the allele |
| `vfdb`, `ecoli_vf`, `victors`, `upec_expec_vf` | `virulence` |
| `plasmidfinder` | `replicon` |
| `bacmet2` | `biocide` |

## Bringing your own database

`gapit db build NAME FASTA` turns any FASTA of reference genes into a fully built
gapit-native database — `records.jsonl`, the `sequences` projection with `gapit/v1`
headers, the BLAST index, and the manifest — in one command.
Here is a two-gene synthetic FASTA plus a metadata TSV (fields below), run against a
scratch datadir:

```console
$ gapit db build tinyamr my_genes.fa --datadir ./db --tsv my_meta.tsv
gapit: generated /tmp/opencode/gapit-build-demo/db/tinyamr/sequences
gapit: self-check passed for tinyamr
gapit: BLAST index built (nucl)
{"db":"tinyamr","records":2,"dbtype":"nucl","destination":"/tmp/opencode/gapit-build-demo/db/tinyamr"}
```

stderr carries the per-step progress, stdout the one-line JSON receipt (same fields
as `db fetch`, see the table above). The database screens immediately:

```console
$ gapit screen --datadir ./db --db tinyamr contig.fa
Processing: contig.fa
Found 1 genes in contig.fa
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
contig.fa	contig1	1	240	+	syn_betalac	1-240/240	===============	0/0	100.00	100.00	tinyamr	SYN-0001	synthetic class A beta-lactamase	ampicillin;cephalosporin
```

The `ACCESSION` and `RESISTANCE` values came from the TSV merge; the `PRODUCT` from
the FASTA description. `--dbtype nucl|prot` forces the molecule type; by default the
abricate mol-type heuristic decides from the sequences themselves. Input may be plain,
`.gz`, or `.bz2`.

### Header detection

The header kind is detected per record, so mixed files work:

| FASTA header | Gene | Accession | Function | Product |
|---|---|---|---|---|
| `>syn_betalac synthetic class A beta-lactamase` | `syn_betalac` | — | — | `synthetic class A beta-lactamase` |
| `>olddb~~~sul1~~~U12338.4:1-940~~~SULFONAMIDE sulfonamide resistance` | `sul1` | `U12338.4:1-940` | `SULFONAMIDE` | `sulfonamide resistance` |
| `>gapit\|db=old\|gene=sul2\|acc=X\|func=streptomycin aminoglycoside` | `sul2` | `X` | `streptomycin` | `aminoglycoside` |

Plain headers without any description text fall back to `--description TEXT`, then
to the gene name. The `db` field of every record is always NAME (the database being
built), and `source_id` keeps the original id token verbatim. Malformed `gapit|`
headers fail the build with `HEADER_MALFORMED` (exit 4); structurally invalid FASTA
fails with `INVALID_FASTA` (exit 5).

### Merging metadata from a TSV

`--tsv FILE` adds or overrides accession and function classes per gene. The rules:

- The header row is mandatory and must contain a `gene` column; `accession` and
  `function` are optional per file (absent columns are simply not merged), and extra
  columns are ignored. A missing `gene` column fails with exit 5
  `METADATA_MALFORMED`.
- Rows are keyed by gene: the **first row wins** on duplicates, with a warning on
  stderr (`--quiet` silences it).
- The `function` column is `;`-separated for multiple classes, e.g.
  `ampicillin;cephalosporin`.
- Genes present only in the TSV produce a stderr warning and are skipped — the FASTA
  is the truth for what exists.

The example TSV used above, in full:

```console
$ cat my_meta.tsv
gene	accession	function
syn_betalac	SYN-0001	ampicillin;cephalosporin
```

### Rebuilding

Like `db fetch`, an existing database is never overwritten silently:

```console
$ gapit db build tinyamr my_genes.fa --datadir ./db
{"schema":"gapit.error/1","code":"DB_ALREADY_EXISTS","message":"won't overwrite existing database tinyamr (use --force)","context":{"db":"tinyamr"}}
$ gapit db build tinyamr my_genes.fa --datadir ./db --tsv my_meta.tsv --force
gapit: generated /tmp/opencode/gapit-build-demo/db/tinyamr/sequences
gapit: self-check passed for tinyamr
gapit: BLAST index built (nucl)
{"db":"tinyamr","records":2,"dbtype":"nucl","destination":"/tmp/opencode/gapit-build-demo/db/tinyamr"}
```

Because `records.jsonl` is the editable truth, editing it (or the FASTA) and
rebuilding with `--force` regenerates every downstream artifact.

### The manual (abricate-style) path

`db build` is the recommended route, but any abricate-format datadir also works
as-is: point `--datadir` at a directory containing `<name>/sequences` (a nucleotide
FASTA with `~~~` headers) and run `gapit setupdb` once to build the BLAST index. The
mol-type heuristic from abricate picks `nucl` vs `prot`.

End to end with a tiny three-gene database copied into a scratch datadir:

```console
$ ls /tmp/opencode/gapit-own-db/db/tinyamr
sequences
$ gapit screen --datadir /tmp/opencode/gapit-own-db/db --db tinyamr /tmp/opencode/gapit-own-db/gap.fa
Processing: /tmp/opencode/gapit-own-db/gap.fa
{"schema":"gapit.error/1","code":"DATABASE_NOT_INDEXED","message":"Database /tmp/opencode/gapit-own-db/db/tinyamr/sequences is not indexed, please try: gapit setupdb","context":{"db":"/tmp/opencode/gapit-own-db/db/tinyamr/sequences"}}
$ gapit setupdb --datadir /tmp/opencode/gapit-own-db/db
Indexed tinyamr (3 sequences, nucl)
$ gapit screen --datadir /tmp/opencode/gapit-own-db/db --db tinyamr /tmp/opencode/gapit-own-db/gap.fa
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
/tmp/opencode/gapit-own-db/gap.fa	contig1	1	97	+	sul1	1-94/94	========/======	1/3	100.00	96.91	tinyamr	U12338.4:1-940	sulfonamide-resistant dihydropteroate synthase Sul1	SULFONAMIDE
```

`gapit setupdb` indexes every subdirectory with a readable `sequences` file and exits 0
when all of them have `sequences.nin` or `sequences.pin`.

### Reads mode and legacy datadirs

For reads screening (see [reads.md](./reads.md)), minimap2 always indexes the `sequences`
FASTA in memory; no `.mmi` is built, and one left in a datadir by an older gapit is never
consulted. Persisted indexes were tried and rejected: a default-built `.mmi` overrides the
`-x sr` preset's indexing parameters (minimap2 warns "-k, -w or -H overridden by prebuilt
index"), which misassigns close homologs (CTX-M/SHV allele divergence in the deciding
benchmark) and ran slower than in-memory indexing at current database scale. Legacy
abricate datadirs take the same in-memory path.

A complete `db build` walkthrough with worked examples for every header format, protein
databases, metadata TSVs, and a troubleshooting table lives in
[Custom databases](./custom-db.md).
