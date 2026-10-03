# Outputs

Every gapit command follows one contract: data on stdout, diagnostics on stderr, errors as
typed JSON envelopes with documented exit codes. This page is the reference for each format.
The schemas are versioned and machine-discoverable, so agents can consume them without
reading this page. See also [screen.md](./screen.md), [reads.md](./reads.md), and
[summary.md](./summary.md) for the commands that produce these outputs.

## stdout and stderr

- **stdout carries data only**: the TSV/CSV table, the JSON document, the Markdown report,
  or the one-line receipts printed by `gapit db fetch` and `gapit db install`. With
  `gapit screen --output PATH` the data goes to the file instead and stdout stays empty.
- **stderr carries diagnostics**: progress lines like `Processing: <file>` during screening,
  fetch progress, warnings. `--quiet` silences stderr; it never touches stdout.
- Output is deterministic: stable sort orders, fixed tool parameters, no wall-clock
  timestamps inside data payloads (only the `created_at` metadata field).
- **Emission timing** (`gapit screen`): tsv/csv/md stream — the header or static frontmatter
  first, then each file's rows/section the moment that file finishes (input order;
  head-of-line under `--jobs`). json is a single document written once at the end. The
  concatenated bytes are identical either way.

## TSV (default format)

`gapit screen` prints abricate-compatible TSV, one header row then one row per surviving
hit. Real output from a three-gene fixture database:

```console
$ gapit screen --datadir /tmp/opencode/gapit-own-db/db --db tinyamr /tmp/opencode/gapit-own-db/gap.fa
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
/tmp/opencode/gapit-own-db/gap.fa	contig1	1	97	+	sul1	1-94/94	========/======	1/3	100.00	96.91	tinyamr	U12338.4:1-940	sulfonamide-resistant dihydropteroate synthase Sul1	SULFONAMIDE
```

### The 15 columns

| # | Column | Type | Meaning |
|---|---|---|---|
| 1 | `FILE` | string | Input path as given on the command line (basename only with `--nopath`) |
| 2 | `SEQUENCE` | string | Contig identifier within the input file |
| 3 | `START` | integer | Alignment start on the contig, 1-based |
| 4 | `END` | integer | Alignment end on the contig, always greater than `START` |
| 5 | `STRAND` | string | `+` or `-`; the gene's strand on the subject, not the read direction |
| 6 | `GENE` | string | Gene name from the database header |
| 7 | `COVERAGE` | string | `sstart-send/slen`: aligned span of the gene over its full length |
| 8 | `COVERAGE_MAP` | string | 15-character alignment sketch, see below |
| 9 | `GAPS` | string | `gapopen/gaps`: gap openings over total gap columns |
| 10 | `%COVERAGE` | float | `100 * (length - gaps) / slen`, displayed with two decimals |
| 11 | `%IDENTITY` | float | The `pident` value printed by BLAST, never post-filtered, two decimals |
| 12 | `DATABASE` | string | Source database name |
| 13 | `ACCESSION` | string | Accession or coordinate span from the database header |
| 14 | `PRODUCT` | string | Product description; `n/a` when the header has none |
| 15 | `RESISTANCE` | string | Resistance or functional category from the database header |

Rows are sorted by `SEQUENCE` (lexicographic) then `START` (numeric) within each input file,
and files appear in argument order.

### The two percentages

Both come straight from the BLAST row, with no smoothing:

- **%IDENTITY** is BLAST's `pident`, displayed `%.2f`. It is never re-filtered: `--minid`
  is enforced inside blastn via `-perc_identity`.
- **%COVERAGE** is `100 * (length - gaps) / slen`, the ungapped aligned columns divided by
  the full gene length. The `--mincov` threshold compares the unrounded float, then the
  value is displayed `%.2f`. A hit at 79.996% displays as `80.00` but is discarded. This is
  abricate's exact behavior, kept on purpose.

### COVERAGE_MAP

The minimap sketches where the alignment lands on the gene in 15 character cells:

- Normally there are 15 cells. When the alignment has gap openings, the last cell is given
  to a `/` marker, leaving 14 cells.
- Each cell covers `slen / cells` gene bases. Cells between the alignment's start and end
  get `=`, everything else `.`. With gaps, `/` sits in the middle cell position.
- Coordinates are truncated to integers, so edge cells can stay `.` even for a full-length
  alignment on a long gene. abricate has the same quirk; gapit reproduces it rather than
  "fixing" it.

The example row above is a gapped alignment of `sul1` (94 bp gene, 1 gap opening spanning 3
columns):

```
========/======
|||||||| ||||||
01234567 89...13   cells 0..13 (14, because gapopen > 0)
                   '/' sits after cell 7, the middle position
```

All 14 cells show `=` because the alignment spans the whole gene; the `/` records that the
alignment is gapped. An ungapped full-length alignment prints 15 `=` characters, like this
`tetA` row:

```text
/tmp/opencode/gapit-own-db/full.fa	contig1	1	79	+	tetA	1-79/79	===============	0/0	100.00	100.00	tinyamr	NC_000913.3:100-900	tetracycline efflux pump TetA	TETRACYCLINE
```

### Variations

- `--format csv` switches the separator to `,`.
- `--noheader` drops the `#FILE ...` header row.
- `--nopath` basenames the `FILE` column.
- Overlapping genes at different query spans are both reported. Hits sharing one
  `(contig, start, end)` deduplicate, first BLAST row wins. No interval merging happens on
  the default path; this matches abricate.

### Typed gene databases: designation is a second command

A gene database carrying a `typing.json` (a **typed** database — the bundled `ecoli_dec` is
the example) screens exactly like an untyped one: the frozen 15-column abricate table above,
byte-for-byte, in every format. Designation is the second stage of the two-stage pipeline —
[`gapit typing`](./typing.md#the-two-stage-designation-workflow) reads a screen result
table and renders the calls (real output; `screen --output` wrote the table):

```console
$ gapit screen dec_s2_pic_astA_uidA.fasta --db ecoli_dec --output dec.tsv --quiet
$ gapit typing dec.tsv --quiet
FILE	SCHEME	PHENOTYPE	GENES	CONFIDENCE	SCORE	RUNNER_UP	NOTES
dec_s2_pic_astA_uidA.fasta	gb4789_6	EAEC	astA;pic;uidA	high	1.0000	EHEC (0.0000)	GB 4789.6: any of aggR/pic/astA
dec_s2_pic_astA_uidA.fasta	risk_monitoring	non-DEC	astA;pic;uidA	low	0.0000	STEC (0.0000)	GB 4789.6-2016: EHEC = stx (any subunit) + escV (LEE); severity order EHEC>STEC/EPEC>ETEC>EIEC>EAEC: rules are declared in severity order so hybrids surface as runner_up
```

The pic+astA profile without aggR is the headline divergence — GB 4789.6 calls EAEC while
the risk-monitoring scheme falls back to non-DEC. Zero-hit files emit no screen rows, so
they cannot appear in a typing run (a table with no data rows at all is the
`TYPING_NO_DATA` input error).

## JSON

The screening surfaces emit versioned JSON documents, all requested with `--format json`.
Every document starts with a `schema` field naming its contract. Field names are snake_case
with explicit units (`identity_pct`, `coverage_pct`, `breadth_pct`). Numeric hit values are
rounded to two decimals, matching what the TSV displays.

### gapit.report/1 (contig screening)

Real document, same run as the TSV above:

```json
{
  "schema": "gapit.report/1",
  "tool": {
    "name": "gapit",
    "version": "0.5.4"
  },
  "created_at": "2026-09-19T01:12:20Z",
  "params": {
    "db": "tinyamr",
    "minid": 80.0,
    "mincov": 80.0,
    "threads": 1
  },
  "files": [
    {
      "file": "/tmp/opencode/gapit-own-db/gap.fa",
      "hits": [
        {
          "sequence": "contig1",
          "start": 1,
          "end": 97,
          "strand": "+",
          "gene": "sul1",
          "coverage": "1-94/94",
          "coverage_map": "========/======",
          "gaps": "1/3",
          "coverage_pct": 100.0,
          "identity_pct": 96.91,
          "database": "tinyamr",
          "accession": "U12338.4:1-940",
          "product": "sulfonamide-resistant dihydropteroate synthase Sul1",
          "resistance": "SULFONAMIDE"
        }
      ]
    }
  ]
}
```

| Field | Type | Meaning |
|---|---|---|
| `schema` | string | Always `gapit.report/1` |
| `tool` | object | `{name, version}` of gapit |
| `created_at` | string | ISO-8601 UTC timestamp, second precision |
| `params` | object | Screening parameters in effect |
| `params.db` | string | Database name |
| `params.minid` | number | Minimum identity threshold |
| `params.mincov` | number | Minimum coverage threshold |
| `params.threads` | integer | BLAST thread count |
| `files` | array | One entry per input file, zero-hit files included |
| `files[].file` | string | Input path as given |
| `files[].hits` | array | Surviving hits, empty list when none |

Each hit mirrors the TSV columns one to one:

| Field | Type | Meaning |
|---|---|---|
| `sequence` | string | Contig identifier |
| `start`, `end` | integer | Contig alignment span, 1-based |
| `strand` | string | `+` or `-` |
| `gene` | string | Gene name |
| `coverage` | string | `sstart-send/slen` |
| `coverage_map` | string | The 15-character minimap |
| `gaps` | string | `gapopen/gaps` |
| `coverage_pct` | number | Percent coverage, two decimals |
| `identity_pct` | number | Percent identity, two decimals |
| `database` | string | Source database name |
| `accession` | string | Accession or coordinate span |
| `product` | string | Product description |
| `resistance` | string | Resistance or functional category (frozen name, kept for TSV parity) |

Two additive optional fields appear on hits produced by `gapit screen --merge-fragments`
(cross-contig fragment merging; [screen.md](./screen.md#fragment-merging-merge-fragments))
and are absent from every other row:

| Field | Type | Meaning |
|---|---|---|
| `merged` | boolean | Always `true`; marks a hit assembled from >= 2 gene fragments |
| `fragments` | array | One `{contig, start, end, strand, identity_pct, coverage_pct}` entry per contributing fragment, sorted by `(contig, start)` |

The merged hit itself carries the union `%COVERAGE`, the aligned-length-weighted mean
`%IDENTITY`, summed `GAPS`, and the anchor fragment's query coordinates; `sequence` is the
comma-joined contig list. These fields are additive to `gapit.report/1` (no rename, no
retype, no schema version bump).

### gapit.typing_result/1 (designation from screen results)

`gapit typing RESULT.tsv [RESULT2.tsv ...]` reads one or more gapit/abricate screen result
tables (TSV or CSV), resolves the database from the rows' shared `DATABASE` column, and
evaluates its `typing.json` over each FILE's genes — every `(FILE, GENE)` folds to its best
row by `(%IDENTITY, %COVERAGE)`. Real document (trimmed):

```json
{
  "schema": "gapit.typing_result/1",
  "tool": {"name": "gapit", "version": "0.5.4"},
  "created_at": "2026-10-02T10:04:55Z",
  "source": ["dec.tsv"],
  "db": "ecoli_dec",
  "files": [
    {
      "file": "dec_s2_pic_astA_uidA.fasta",
      "genes": ["astA", "pic", "uidA"],
      "phenotypes": {
        "gb4789_6": {
          "phenotype": "EAEC",
          "score": 1.0,
          "confidence": "high",
          "components": [
            {"name": "aggR", "score": 0.0},
            {"name": "pic", "score": 1.0},
            {"name": "astA", "score": 1.0}
          ],
          "runner_up": {"phenotype": "EHEC", "score": 0.0},
          "notes": ["GB 4789.6: any of aggR/pic/astA"]
        },
        "risk_monitoring": {
          "phenotype": "non-DEC",
          "score": 0.0,
          "confidence": "low",
          "components": [
            {"name": "escV", "score": 0.0},
            {"name": "stx1a", "score": 0.0},
            {"name": "stx1b", "score": 0.0},
            {"name": "stx2a", "score": 0.0},
            {"name": "stx2b", "score": 0.0}
          ],
          "runner_up": {"phenotype": "STEC", "score": 0.0},
          "notes": [
            "GB 4789.6-2016: EHEC = stx (any subunit) + escV (LEE)",
            "severity order EHEC>STEC/EPEC>ETEC>EIEC>EAEC: rules are declared in severity order so hybrids surface as runner_up"
          ]
        }
      }
    }
  ]
}
```

| Field | Type | Meaning |
|---|---|---|
| `schema` | string | Always `gapit.typing_result/1` |
| `tool` | object | `{name, version}` of gapit |
| `created_at` | string | ISO-8601 UTC timestamp, second precision |
| `source` | array | The result table path(s) as given |
| `db` | string | The database every row screened against |
| `files` | array | One entry per FILE that has data rows, first-appearance order |
| `files[].genes` | array | That FILE's present gene names, sorted — the hits behind the calls |
| `files[].phenotypes` | object | Scheme name → its call (the additive shape below) |

Each call carries `phenotype` (null on an ambiguous call), `score`, `confidence`
(`high`/`ambiguous`/`low`), `components[{name, score}]`, and the optional `runner_up`,
`ambiguous[]` (the tied pair), and `notes[]` — the same bodies the evaluator produces on
the cluster path. The TSV/Markdown projection flattens each call to the eight `FILE`,
`SCHEME`, `PHENOTYPE`, `GENES`, `CONFIDENCE`, `SCORE`, `RUNNER_UP`, `NOTES` columns:
ambiguous calls render `-` for the phenotype and carry the candidate pair in `NOTES` (the
runner-up cell also renders `-` there — the pair already speaks); `GENES` follows the
phenotype with the FILE's sorted `;`-joined gene list, repeated on every scheme row of that FILE.
Typed errors: mixed `DATABASE` values across rows are `DATABASE_MISMATCH`, a table with no
data rows `TYPING_NO_DATA`,
a database without a `typing.json` `TYPING_NO_SCHEME`, and a cluster database
`TYPING_CLUSTER_DB` (its typing is integrated into `gapit screen`). Its JSON Schema prints
with `gapit schema typing_result`.

### gapit.reads/1 (reads and assembly screening)

Reads mode reports per-gene coverage across the read set (or assembly; `--r1` accepts FASTQ and
FASTA) instead of per-hit rows. Gene entries sort by `breadth_pct` descending; genes with zero
mapped reads are omitted. Fields from `gapit schema reads`:

| Field | Type | Meaning |
|---|---|---|
| `schema` | string | Always `gapit.reads/1` |
| `tool` | object | `{name, version}` of gapit |
| `created_at` | string | ISO-8601 UTC timestamp, second precision |
| `params` | object | Read-screening parameters in effect |
| `params.db` | string | Database name |
| `params.read_type` | string | minimap2 preset in effect: `sr`, `map-ont`, or `map-hifi` (resolved from the detected input) |
| `params.min_breadth` | number | Presence threshold, default 90.0 |
| `params.threads` | integer | minimap2 thread count |
| `files` | array | One entry per read set |
| `files[].reads` | array of string | Input read/assembly paths (all lanes) |
| `files[].genes` | array | Per-gene presence calls |

Each gene entry:

| Field | Type | Meaning |
|---|---|---|
| `gene` | string | Gene name |
| `database` | string | Source database name |
| `accession` | string | Accession or coordinate span |
| `product` | string | Product description |
| `resistance` | string | Functional category |
| `tlen` | integer | Gene length in bases |
| `breadth_pct` | number | Percent of gene bases covered by at least one primary alignment |
| `mean_depth` | number | Mean per-base depth over the gene length |
| `reads_mapped` | integer | Distinct reads with a primary alignment on the gene |
| `present` | boolean | `breadth_pct >= min_breadth` |

### gapit.reads/2 (reads screening with identity/MAPQ filtering)

Opt-in variant of reads/1, emitted only when `--min-identity` or `--min-mapq` is nonzero: PAF
alignments below the thresholds are dropped before aggregation, which removes the
family-splitting over-calls breadth-only presence suffers on homologous genes. With both
thresholds off the output stays `gapit.reads/1`, byte-identical. Fields from
`gapit schema reads2` — identical to reads/1 except:

| Field | Type | Meaning |
|---|---|---|
| `schema` | string | Always `gapit.reads/2` |
| `params.min_identity` | number | Per-alignment identity floor in effect, 0 = off |
| `params.min_mapq` | integer | Per-alignment MAPQ floor in effect, 0 = off |
| `files[].genes[].mean_identity_pct` | number | Alignment-length-weighted mean per-alignment identity over the kept alignments, rounded to 2 |

Per-alignment identity is `100 * (alen - nm) / alen` (PAF block length and `NM:i:` tag; a row
without NM counts as 100). See [reads.md](./reads.md#filtering-alignments-by-identity-and-mapq-gapitreads2)
for the homolog worked example and threshold guidance.

### gapit.cluster/1 (cluster-database screening)

Screening a `kind: cluster` database (minimap2 `asm20` engine; [screen.md](./screen.md#cluster-databases-kind-cluster))
emits one document per run: the best-locus call per file, every covered locus with its
per-gene verdicts, and — on databases carrying a `typing.json` — the phenotype call with its
explainable score breakdown. Fields from `gapit schema cluster`:

| Field | Type | Meaning |
|---|---|---|
| `schema` | string | Always `gapit.cluster/1` |
| `tool` | object | `{name, version}` of gapit |
| `created_at` | string | ISO-8601 UTC timestamp, second precision |
| `params` | object | `db`, `preset` (always `asm20`), `min_gene_cov`, `min_gene_id`, `min_cluster_cov`, `threads` |
| `files` | array | One entry per input file |
| `files[].file` | string | Input path as given |
| `files[].best` | object \| null | The best-locus call, null when no locus cleared `min_cluster_cov` |
| `files[].loci` | array | Every covered locus in rank order (coverage desc, identity desc, covered bp desc, id asc) |

The `best` block:

| Field | Type | Meaning |
|---|---|---|
| `best.locus` | string | Rank-1 locus id |
| `best.label` / `best.type` | string | The locus label and type from the source notes (kaptive `K locus`/`K type` semantics) |
| `best.coverage_pct` / `best.identity_pct` | number | Union coverage and matches-weighted identity over all records, two decimals |
| `best.genes_present` / `genes_partial` / `genes_absent` | integer | Verdict counts for the best locus's annotated genes |
| `best.phenotype` | string \| null | The typing call; null when ambiguous or on untyped databases |
| `best.phenotype_detail` | object | Additive, present only on typed runs: `{score, confidence, components[{name, score}], runner_up, ambiguous[]}` |

Each `loci[]` entry carries `locus`, `label`, `type`, `coverage_pct`, `identity_pct`,
`rank`, `missing` (gene ids not fully present), and `genes[]` with per-gene `start`/`end`
(1-based, locus coordinates), `strand`, `coverage_pct`, `identity_pct`, and `verdict`
(`present`/`partial`/`absent`).

`phenotype_detail.confidence` is `high` (winner ≥ cutoff and separated by ≥
`ambiguity_margin`), `ambiguous` (two rules inside the margin: phenotype null, the top two
listed in `ambiguous`), or `low` (below cutoff: the document's fallback string). Component
names are the scored rule's own: gene ids for `weighted_genes`, `coverage`/`identity`/
`key_genes` for `cluster_match`, feature strings (plus `bias`) for `learned_linear`.

The TSV form is one row per file with the header `FILE BEST_LOCUS TYPE PHENOTYPE COVERAGE
IDENTITY PRESENT PARTIAL MISSING_IDS` on typed databases (`PHENOTYPE` is `-` when
ambiguous or uncalled) and the same header without `PHENOTYPE` on untyped ones. A file with
no locus call renders `- - - 0.00 0.00 0 0 -`-shaped dashes.

### gapit.features/1 and gapit.typing/1 (database-side documents)

Two cluster-database artifacts are versioned documents you can introspect with `gapit
schema` but that never appear on stdout: `features` (the `features.json` locus/gene feature
table written by `db build`) and `typing` (the declarative scoring spec installed with
`db build --typing FILE` — rules of kind `weighted_genes` / `cluster_match` /
`learned_linear`, plus `cutoff`, `ambiguity_margin`, and `fallback`). The typing document is
the contract the phenotype evaluator implements; see [databases.md](./databases.md#cluster-databases-gbkgff).

### gapit.summary/1 (summary matrix)

`gapit summary --format json` turns report tables into a gene matrix. Real output
summarizing two report files, one with a `tetA` hit and one with none:

```json
{
  "schema": "gapit.summary/1",
  "tool": {
    "name": "gapit",
    "version": "0.5.4"
  },
  "created_at": "2026-09-19T01:11:09Z",
  "params": {
    "metric": "%COVERAGE",
    "nopath": false
  },
  "genes": [
    "tetA"
  ],
  "rows": [
    {
      "file": "/tmp/opencode/gapit-docs-Ifvc44/full.tsv",
      "num_found": 1,
      "cells": {
        "tetA": [
          "100.00"
        ]
      }
    },
    {
      "file": "/tmp/opencode/gapit-docs-Ifvc44/partial.tsv",
      "num_found": 0,
      "cells": {}
    }
  ]
}
```

| Field | Type | Meaning |
|---|---|---|
| `schema` | string | Always `gapit.summary/1` |
| `tool` | object | `{name, version}` of gapit |
| `created_at` | string | ISO-8601 UTC timestamp, second precision |
| `params.metric` | string | `%COVERAGE` or `%IDENTITY`, per the `--identity` flag |
| `params.nopath` | boolean | Whether row labels were basenamed |
| `genes` | array of string | Sorted union of all gene names across inputs |
| `rows` | array | One row per input file |
| `rows[].file` | string | Input report path (the row label) |
| `rows[].num_found` | integer | Count of distinct genes found in that report |
| `rows[].cells` | object | Gene name to list of original cell strings; absent genes omitted |

Cell values keep the exact strings from the input reports, in file order, so `%COVERAGE`
and `%IDENTITY` summaries round-trip without reformatting.

### Stability policy

Schema names are semver'd (`gapit.report/1`). Within a major version, existing fields are
never renamed or retyped; additive changes come with a schema-version bump and a note in
`PLAN.md`. Validate anything you build against the schemas printed by `gapit schema`.

## Markdown

`--format md` renders the same data as a human-readable Markdown report: a STATIC YAML
frontmatter (`schema`, `tool`, `created_at`, `db`, `minid`, `mincov`, `threads` — no run
totals, which live in the JSON document), then one section per input file containing the
same 15 columns as the TSV in a pipe table. The frontmatter is knowable before the first
file, so md streams per file exactly like the tsv rows. Reads mode and summary mode have
analogous Markdown forms. Frontmatter gives parsers a stable header; the tables read
naturally in a terminal or editor.

## Errors

Any failure prints one line of JSON to stderr:

```text
{"schema": "gapit.error/1", "code": "<CODE>", "message": "...", "context": {...}}
```

| Field | Type | Meaning |
|---|---|---|
| `schema` | string | Always `gapit.error/1` |
| `code` | string | Machine-stable error code, snake_case |
| `message` | string | Human-readable description |
| `context` | object | String-valued details (file paths, digests, db names) |

The exit code tells you the failure class:

| Exit code | Meaning |
|---|---|
| 0 | success |
| 1 | unexpected error (envelope code `UNEXPECTED`) |
| 2 | usage error |
| 3 | missing dependency |
| 4 | database error |
| 5 | input error |

Screening a file that doesn't exist exits 5 with this real stderr line:

```console
$ gapit screen /nonexistent/contigs.fa --db ncbi
{"schema":"gapit.error/1","code":"INPUT_NOT_FOUND","message":"input file not found or unreadable: /nonexistent/contigs.fa","context":{"file":"/nonexistent/contigs.fa"}}
$ echo $?
5
```

Two more, each produced by the command shown:

```console
$ gapit db fetch nosuchdb            # exit 2, usage
{"schema":"gapit.error/1","code":"USAGE_ERROR","message":"unknown database: nosuchdb (available: argannot, bacmet2, card, ecoh, ecoli_vf, megares, ncbi, plasmidfinder, resfinder, upec_expec_vf, vfdb, victors)","context":{"db":"nosuchdb"}}

$ gapit screen --datadir /tmp/opencode/gapit-own-db/db --db tinyamr contigs.fa   # exit 4, db error, before indexing
{"schema":"gapit.error/1","code":"DATABASE_NOT_INDEXED","message":"Database /tmp/opencode/gapit-own-db/db/tinyamr/sequences is not indexed, please try: gapit setupdb","context":{"db":"/tmp/opencode/gapit-own-db/db/tinyamr/sequences"}}
```

A missing external binary (exit 3) reports which one:

```console
{"schema":"gapit.error/1","code":"MISSING_DEPENDENCY","message":"required binary not found on PATH: blastdbcmd","context":{"binary":"blastdbcmd"}}
```

stdout stays empty on failure: no partial tables, no partial JSON.

## Self-description

An agent can discover the whole contract from the binary alone.

```console
$ gapit --version --json
{"schema":"gapit.version/1","name":"gapit","version":"0.5.4"}
```

`gapit schema <name>` prints the JSON Schema for each document. The eleven names: `report`,
`typing_result`, `reads`, `reads2`, `cluster`, `summary`, `error`, `version`, plus the
database-side documents `features` (gapit.features/1, a cluster db's feature table),
`typing` (gapit.typing/1, a cluster db's declarative scoring spec), and `floors`
(gapit.floors/1, per-gene identity floors for reads-mode presence). A trimmed fragment of
`gapit schema report`:

```json
{
  "description": "gapit.report/1 \u2014 the canonical machine-readable screening output.",
  "properties": {
    "schema": {
      "const": "gapit.report/1",
      "default": "gapit.report/1",
      "title": "Schema",
      "type": "string"
    },
    "created_at": {
      "title": "Created At",
      "type": "string"
    },
    "params": {
      "$ref": "#/$defs/ParamsDocument"
    },
    "files": {
      "items": {
        "$ref": "#/$defs/FileDocument"
      },
      "title": "Files",
      "type": "array"
    }
  },
  "required": [
    "created_at",
    "params",
    "files"
  ],
  "title": "ReportDocument",
  "type": "object"
}
```

Point any JSON Schema validator at this output to check a document before consuming it.
The [agent guide](./agents.md) shows the full introspection workflow.
