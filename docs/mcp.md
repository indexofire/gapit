# MCP server

gapit ships a Model Context Protocol (MCP) stdio server that exposes the CLI as
tools. An agent runtime can screen assemblies, summarize reports, introspect
schemas, and self-provision databases without shelling out or parsing
terminal output. The protocol is hand-rolled newline-delimited JSON-RPC 2.0
with no extra dependencies: one request per line on stdin, one response per
line on stdout.

## Entry points

| Command | What it is |
|---|---|
| `gapit mcp` | The `mcp` subcommand of the CLI. |
| `gapit-mcp` | Console script installed with the package (same code). |

Both serve until stdin closes (EOF ends the process cleanly). stdout carries
protocol lines only; stderr is reserved for protocol-internal errors.

## Client configuration

Register the console script with an MCP client, Claude-style:

```json
{"mcpServers": {"gapit": {"command": "gapit-mcp"}}}
```

If `gapit-mcp` is not on the client's PATH, use the absolute path, or invoke
the module CLI as `{"command": "gapit", "args": ["mcp"]}`.

## Tools

Eight tools: four read-only analysis tools and four database
self-provisioning tools (`db_fetch`, `db_build`, `db_search`,
`db_outdated`) that mirror `gapit db ...` exactly. Argument lists below are
transcribed from the live `tools/list` `inputSchema` objects.

| Tool | Arguments | Returns |
|---|---|---|
| `screen` | `files` (array of strings, required), `db` (string, default `ncbi`), `minid` (number), `mincov` (number), `format` (string: `json` \| `tsv` \| `md`, default `json`) | Report text: `gapit.report/1` JSON by default, TSV or Markdown per `format` |
| `summary` | `files` (array of report table paths, required), `identity` (boolean), `nopath` (boolean) | `gapit.summary/1` JSON |
| `schema` | `name` (string, required, one of `error`, `list`, `reads`, `report`, `summary`, `version`) | The JSON Schema of that output document |
| `db_list` | none | `gapit.dblist/1`: provider names, install state, record counts |
| `db_fetch` | `name` (string), `datadir` (string), `force` (boolean, default `false`) | One JSON receipt line per database (`db`, `records`, `dbtype`, `destination`). Name omitted: the card+vfdb default set from bundled snapshots. Network installs can take minutes |
| `db_build` | `name` (string, required), `fasta` (string, required — a LOCAL filesystem path), `tsv` (string), `dbtype` (string: `nucl` \| `prot`), `description` (string), `datadir` (string), `force` (boolean, default `false`) | One JSON receipt line (`db`, `records`, `dbtype`, `destination`) |
| `db_search` | `term` (string, required), `db` (string), `field` (string: `gene` \| `accession` \| `function` \| `product` \| `any`, default `any`), `exact` (boolean, default `false`), `limit` (integer ≥ 0, default `100`; `0` = unlimited), `datadir` (string) | TSV hit rows with columns `DB`, `GENE`, `ACCESSION`, `FUNCTION`, `PRODUCT`, `LENGTH` |
| `db_outdated` | `days` (integer ≥ 0, default `90`), `datadir` (string) | TSV rows with columns `NAME`, `FETCHED_AT`, `AGE_DAYS`, `STATUS` |

Omitted `screen` thresholds fall back to the CLI defaults (`minid` 80,
`mincov` 80), matching `gapit screen`.

Notes:

- File paths resolve relative to the server's working directory; absolute
  paths are safest.
- The datadir comes from the `GAPIT_DATADIR` environment variable, then
  `~/.local/share/gapit/db` (see `./databases.md`). The db tools also accept
  an explicit `datadir` argument per call.
- `screen` in MCP covers contigs only; reads mode and `csv` are CLI-only.
- `db_build`'s `fasta` must be a path on the server's filesystem — the agent
  provides a local path, not file contents.
- Tool failures do not use JSON-RPC errors. They return `isError: true` with
  the `gapit.error/1` envelope serialized as the text content.

## Example session

Build a throwaway datadir from the repo's test fixture (same recipe the test
suite uses), then talk to the server directly:

```bash
export PATH="$PWD/.pixi/envs/default/bin:$PATH"
export GAPIT_DATADIR=/tmp/gapit-mcp-demo/datadir
mkdir -p "$GAPIT_DATADIR"
cp -r tests/data/db/tinyamr "$GAPIT_DATADIR/"
makeblastdb -in "$GAPIT_DATADIR/tinyamr/sequences" -dbtype nucl \
  -out "$GAPIT_DATADIR/tinyamr/sequences" > /dev/null
```

Feed three request lines to `gapit mcp`:

```bash
cat <<'EOF' | gapit mcp
{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18"}}
{"jsonrpc":"2.0","id":2,"method":"tools/list"}
{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"screen","arguments":{"files":["tests/data/contigs/full.fa"],"db":"tinyamr"}}}
EOF
```

Response line 1, verbatim:

```text
{"jsonrpc":"2.0","id":1,"result":{"protocolVersion":"2025-06-18","capabilities":{"tools":{}},"serverInfo":{"name":"gapit","version":"0.1.0"}}}
```

Response line 2 (real output, elided in the middle; each tool carries its full
`inputSchema`):

```text
{"jsonrpc":"2.0","id":2,"result":{"tools":[{"name":"screen","description":"Screen contig files for AMR/virulence genes (json = gapit.report/1).","inputSchema":{"type":"object","properties":{"files":{"type":"array","items":{"type":"string"}},"db":{"type":"string","default":"ncbi"}, ... },"required":["files"]}}, {"name":"summary", ...}, {"name":"schema", ...}, {"name":"db_list", ...}, {"name":"db_fetch", ...}, {"name":"db_build", ...}, {"name":"db_search", ...}, {"name":"db_outdated", ...}]}}
```

Response line 3 (real output, text content elided). The screen result rides in
`result.content[0].text` as a JSON string:

```text
{"jsonrpc":"2.0","id":3,"result":{"content":[{"type":"text","text":"{\n  \"schema\": \"gapit.report/1\", ... }"}],"isError":false}}
```

Decoded, that text is a complete `gapit.report/1` document (real output from
the same run):

```json
{
  "schema": "gapit.report/1",
  "tool": {"name": "gapit", "version": "0.1.0"},
  "created_at": "2026-09-19T01:15:44Z",
  "params": {"db": "tinyamr", "minid": 80.0, "mincov": 80.0, "threads": 1},
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

### Failure shapes

A missing input file returns a tool error, not a protocol error (real lines):

```text
{"jsonrpc":"2.0","id":5,"result":{"content":[{"type":"text","text":"{\"schema\":\"gapit.error/1\",\"code\":\"INPUT_NOT_FOUND\",\"message\":\"input file not found or unreadable: /tmp/gapit-mcp-demo/nope.fa\",\"context\":{\"file\":\"/tmp/gapit-mcp-demo/nope.fa\"}}"}],"isError":true}}
```

Parse `content[0].text` as JSON and branch on `code` (codes and exit-code
mapping: `./outputs.md`).

### Self-provisioning: build a database, then screen

The db tools let an agent provision its own databases mid-session. Real
session against a fresh datadir (`my_genes.fa` holds one 240 bp synthetic
gene, `query.fa` a 200 bp substring of it):

```bash
mkdir -p /tmp/gapit-mcp-demo/datadir
export GAPIT_DATADIR=/tmp/gapit-mcp-demo/datadir
```

Two request lines and their verbatim responses:

```text
{"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"db_build","arguments":{"name":"myamr","fasta":"/tmp/gapit-mcp-demo/my_genes.fa"}}}
{"jsonrpc":"2.0","id":2,"result":{"content":[{"type":"text","text":"{\"db\":\"myamr\",\"records\":1,\"dbtype\":\"nucl\",\"destination\":\"/tmp/gapit-mcp-demo/datadir/myamr\"}"}],"isError":false}}
{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"screen","arguments":{"files":["/tmp/gapit-mcp-demo/query.fa"],"db":"myamr"}}}
{"jsonrpc":"2.0","id":3,"result":{"content":[{"type":"text","text":"{ ... gapit.report/1 ... }"}],"isError":false}}
```

Decoded, the screen text is a complete `gapit.report/1` document whose only
hit is the gene the agent just built the database from (real output):

```json
{
  "schema": "gapit.report/1",
  "tool": {"name": "gapit", "version": "0.1.0"},
  "created_at": "2026-09-20T13:38:39Z",
  "params": {"db": "myamr", "minid": 80.0, "mincov": 80.0, "threads": 1},
  "files": [
    {
      "file": "/tmp/gapit-mcp-demo/query.fa",
      "hits": [
        {
          "sequence": "contig1",
          "start": 1,
          "end": 200,
          "strand": "+",
          "gene": "demov2",
          "coverage": "1-200/240",
          "coverage_map": "=============..",
          "gaps": "0/0",
          "coverage_pct": 83.33,
          "identity_pct": 100.0,
          "database": "myamr",
          "accession": "",
          "product": "demo beta-lactamase variant 2",
          "resistance": ""
        }
      ]
    }
  ]
}
```

The same session pattern works with `db_fetch` (installs card+vfdb from the
bundled snapshots when `name` is omitted) and `db_search`/`db_outdated` for
inspection. Custom database construction rules (header kinds, `--tsv`
metadata): `./custom-db.md`.

## Protocol notes

- **Framing.** One JSON-RPC 2.0 message per stdin line, one response line per
  request on stdout, UTF-8, LF endings. Responses use compact separators.
- **initialize.** Echoes a requested string `protocolVersion`; when absent or
  empty it answers with the built-in default `2025-06-18`. The result carries
  `capabilities.tools` and `serverInfo` (`name: gapit`, package version).
- **Methods.** Exactly three are handled: `initialize`, `tools/list`,
  `tools/call`. Anything else gets a JSON-RPC error `-32601`:

  ```text
  {"jsonrpc":"2.0","id":4,"error":{"code":-32601,"message":"method not found: resources/list"}}
  ```

- **Unknown tool.** `tools/call` naming an unregistered tool, or with invalid
  params, gets `-32602`:

  ```text
  {"jsonrpc":"2.0","id":6,"error":{"code":-32602,"message":"unknown tool: nope"}}
  ```

- **Batch arrays are unsupported.** A JSON array line is not a valid frame and
  is dropped; send one request per line.
- **Notifications are silent.** Frames without an `id` (for example
  `notifications/initialized`) never produce a response.
- **Junk tolerance.** Non-JSON lines are ignored without a parse-error reply,
  and the loop keeps serving. EOF terminates cleanly.
