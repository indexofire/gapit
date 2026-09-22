# gapit for agents

gapit is built to be driven by programs. Every output is machine-readable by
design, the whole contract is discoverable from the binary itself, and failures
are typed JSON instead of prose. This page is the contract; you should not need
any other documentation to integrate.

## Discovery without docs

Three commands answer "what can this thing do" without reading a manual.

Version, as one JSON line:

```console
$ gapit --version --json
{"schema":"gapit.version/1","name":"gapit","version":"0.2.0"}
```

Output schemas. Six documents are introspectable: `report`, `reads`,
`summary`, `list`, `error`, `version`. Each prints its full JSON Schema:

```console
$ gapit schema report      # gapit.report/1 (contig screening)
$ gapit schema reads       # gapit.reads/1  (FASTQ screening)
$ gapit schema summary     # gapit.summary/1
$ gapit schema list        # gapit.list/1
$ gapit schema error       # gapit.error/1
$ gapit schema version     # gapit.version/1
```

For example, `gapit schema version` (real output):

```json
{
  "description": "gapit.version/1 \u2014 compact one-line self-description.",
  "properties": {
    "schema": {
      "const": "gapit.version/1",
      "default": "gapit.version/1",
      "title": "Schema",
      "type": "string"
    },
    "name": {
      "default": "gapit",
      "title": "Name",
      "type": "string"
    },
    "version": {
      "title": "Version",
      "type": "string"
    }
  },
  "required": ["version"],
  "title": "VersionDocument",
  "type": "object"
}
```

Installed databases (real output against a datadir holding one database):

```console
$ gapit list --json
{
  "schema": "gapit.list/1",
  "databases": [
    {
      "name": "tinyamr",
      "sequences": 3,
      "dbtype": "nucl",
      "date": "2026-Sep-19"
    }
  ]
}
```

## Deterministic outputs

Same input, same bytes. That holds because:

- Sort orders are fixed and documented (`./screen.md`, `./summary.md`).
- No wall-clock timestamps inside data payloads. The only time field is
  document metadata (`created_at`).
- LF line endings, UTF-8, stable key order per schema.
- **stdout purity**: data on stdout, diagnostics on stderr, always. Redirect
  stdout to a file and you get exactly the document, nothing else.

## Exit codes as control flow

The integer on exit tells you the failure class before you read a byte of
stderr.

| Exit code | Meaning |
|---|---|
| 0 | success |
| 1 | unexpected error |
| 2 | usage error |
| 3 | missing dependency |
| 4 | database error |
| 5 | input error |

Envelope details and error codes: `./outputs.md`.

## JSON stability policy

- Every document self-identifies with a version string: `gapit.report/1`,
  `gapit.reads/1`, `gapit.summary/1`, `gapit.list/1`, `gapit.error/1`,
  `gapit.version/1`. Check `schema` first, dispatch on it.
- Schemas follow semver. A minor bump never renames or retypes an existing
  field; new fields may appear, so ignore unknown keys rather than rejecting
  them.
- Keys are snake_case with explicit units (`identity_pct`, `coverage_pct`).
- The authoritative machine contract is the schema itself: `gapit schema
  <doc>` is always current, this page may lag.

## Error envelope parsing

Any failure prints exactly one JSON line to stderr:

```json
{"schema":"gapit.error/1","code":"...","message":"...","context":{...}}
```

Real example, screening against a database that does not exist (exit 4):

```text
{"schema":"gapit.error/1","code":"DATABASE_NOT_FOUND","message":"Database nosuchdb is not in /tmp/gapit-mcp-demo/datadir. Available: tinyamr","context":{"db":"nosuchdb","datadir":"/tmp/gapit-mcp-demo/datadir"}}
```

And a missing external binary (exit 3):

```text
{"schema":"gapit.error/1","code":"MISSING_DEPENDENCY","message":"required binary not found on PATH: blastn","context":{"binary":"blastn"}}
```

Parsing guidance:

- Progress diagnostics (for example `Processing: sampleA.fa`) also go to
  stderr. The envelope is the one stderr line that parses as JSON with
  `"schema":"gapit.error/1"`.
- `code` is a stable string enum (`USAGE_ERROR`, `MISSING_DEPENDENCY`,
  `DATABASE_NOT_FOUND`, `INPUT_NOT_FOUND`, ...). Branch on `code`, not on
  `message`.
- `context` carries machine-readable details (file paths, db name), safe to
  log or surface.

## Recommended workflow

Discover, screen, parse, summarize. Session below run against the repo's test
fixture datadir (setup recipe in `./mcp.md`).

1. Discover available databases (`gapit list --json`, or the `db_list` MCP
   tool), pick a `db` name.

2. Screen each sample as JSON and parse hits straight out of the document:

   ```console
   $ gapit screen sampleA.fa --db tinyamr --format json > sampleA.json
   $ gapit screen sampleB.fa --db tinyamr --format json > sampleB.json
   ```

   ```python
   import json

   for path in ("sampleA.json", "sampleB.json"):
       doc = json.load(open(path))
       assert doc["schema"] == "gapit.report/1"
       for f in doc["files"]:
           for hit in f["hits"]:
               print(f["file"], hit["gene"], hit["identity_pct"], hit["coverage_pct"])
   ```

   Real output:

   ```text
   sampleA.fa tetA 100.0 100.0
   ```

3. For a cross-sample matrix, keep the default TSV reports and feed them to
   `gapit summary`:

   ```console
   $ gapit screen sampleA.fa --db tinyamr > sampleA.tsv
   $ gapit screen sampleB.fa --db tinyamr > sampleB.tsv
   $ gapit summary sampleA.tsv sampleB.tsv
   #FILE	NUM_FOUND	tetA
   sampleA.tsv	1	100.00
   sampleB.tsv	0	.
   ```

Screening options: `./screen.md`. Summary semantics: `./summary.md`.

## Tool-runtime integration

If your runtime speaks MCP, register gapit as a server and call `screen`,
`screen_reads`, `summary`, `schema`, and `db_list` as tools instead of
subprocesses: `./mcp.md`. The server also exposes the `db` commands
(`db_fetch`, `db_build`, `db_search`, `db_outdated`), so an agent can
provision the database it needs and screen against it in one session.
