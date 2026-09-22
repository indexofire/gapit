# Outputs

Every gapit command follows one contract: data on stdout, diagnostics on stderr, errors as
typed JSON envelopes with documented exit codes. This page is the reference for each format.
The schemas are versioned and machine-discoverable, so agents can consume them without
reading this page. See also [screen.md](./screen.md), [reads.md](./reads.md), and
[summary.md](./summary.md) for the commands that produce these outputs.

## stdout and stderr

- **stdout carries data only**: the TSV/CSV table, the JSON document, the Markdown report,
  or the one-line receipts printed by `gapit db fetch` and `gapit db install`.
- **stderr carries diagnostics**: progress lines like `Processing: <file>` during screening,
  fetch progress, warnings. `--quiet` silences stderr; it never touches stdout.
- Output is deterministic: stable sort orders, fixed tool parameters, no wall-clock
  timestamps inside data payloads (only the `created_at` metadata field).

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

## JSON

Three screening surfaces emit versioned JSON documents, all requested with `--format json`.
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
    "version": "0.2.1"
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

### gapit.summary/1 (summary matrix)

`gapit summary --format json` turns report tables into a gene matrix. Real output
summarizing two report files, one with a `tetA` hit and one with none:

```json
{
  "schema": "gapit.summary/1",
  "tool": {
    "name": "gapit",
    "version": "0.2.1"
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

`--format md` renders the same data as a human-readable Markdown report: YAML frontmatter
(`schema`, `tool`, `created_at`, `db`, `minid`, `mincov`, `threads`, `files`, `hits`),
then one section per input file containing the same 15 columns as the TSV in a pipe table.
Reads mode and summary mode have analogous Markdown forms. Frontmatter gives parsers a
stable header; the tables read naturally in a terminal or editor.

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
$ gapit screen /nonexistent/contigs.fa
{"schema":"gapit.error/1","code":"INPUT_NOT_FOUND","message":"input file not found or unreadable: /nonexistent/contigs.fa","context":{"file":"/nonexistent/contigs.fa"}}
$ echo $?
5
```

Two more, each produced by the command shown:

```console
$ gapit db fetch nosuchdb            # exit 2, usage
{"schema":"gapit.error/1","code":"USAGE_ERROR","message":"unknown provider: nosuchdb (available: argannot, bacmet2, card, ecoh, ecoli_vf, megares, ncbi, plasmidfinder, resfinder, upec_expec_vf, vfdb, victors)","context":{"provider":"nosuchdb"}}

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
{"schema":"gapit.version/1","name":"gapit","version":"0.2.1"}
```

`gapit list --json` describes installed databases (`gapit.list/1`, first two of twelve
shown, trimmed):

```json
{
  "schema": "gapit.list/1",
  "databases": [
    {
      "name": "argannot",
      "sequences": 2224,
      "dbtype": "nucl",
      "date": "2026-Sep-18"
    },
    {
      "name": "bacmet2",
      "sequences": 746,
      "dbtype": "prot",
      "date": "2026-Sep-18"
    }
  ]
}
```

| Field | Type | Meaning |
|---|---|---|
| `schema` | string | Always `gapit.list/1` |
| `databases[].name` | string | Database name, usable as `--db` |
| `databases[].sequences` | integer | Sequence count in the database |
| `databases[].dbtype` | string | `nucl` or `prot` |
| `databases[].date` | string | Build date, abricate `%d-%b-%Y` format |

`gapit schema <name>` prints the JSON Schema for each document. The six names: `report`,
`reads`, `summary`, `list`, `error`, `version`. A trimmed fragment of `gapit schema report`:

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
