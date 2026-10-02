# gapit documentation

gapit screens contig files and FASTQ reads for known genes — AMR, virulence, serotype,
plasmid replicons, species targets, or any custom database. It is an agent-first Python
reimplementation of [abricate](https://github.com/tseemann/abricate): abricate-compatible
TSV on stdout, plus first-class JSON and Markdown outputs with versioned schemas.

## Pages

| Page | Contents |
|---|---|
| [Installation](./installation.md) | PyPI wheel and pixi paths, external binaries, verification, database bootstrap, shell completions |
| [Quickstart](./quickstart.md) | A complete first session on the bundled test fixture, offline |
| [Screening contigs](./screen.md) | `gapit screen` on FASTA/GBK/EMBL inputs: thresholds, filters, formats |
| [Screening reads](./reads.md) | FASTQ and assembly FASTA through minimap2: `--r1`/`--r2`, presets, breadth-based presence, two-stage survey |
| [Summary](./summary.md) | `gapit summary`: gene presence/absence matrix across report tables |
| [Databases](./databases.md) | The database catalog, `gapit db fetch/list/install`, datadirs, native `gapit/v1` format, cluster databases (GBK/GFF, typing.json, kaptive) |
| [Custom databases](./custom-db.md) | `gapit db build` walkthrough: any FASTA to a screenable database, worked examples |
| [Typing schemes](./typing.md) | The two-stage designation pipeline (`gapit screen` then `gapit typing`) and `gapit.typing/2`: rule/scheme primitives, the six-scheme cookbook (Doumith, ShigaTyper, meningotype, VP O/K, cholerae, DEC), allele probe tricks, calibration |
| [Outputs](./outputs.md) | TSV/JSON/Markdown formats, schemas, error envelopes, exit codes |
| [MCP server](./mcp.md) | `gapit mcp`: analysis + database tools for agent runtimes over stdio JSON-RPC |
| [Agent guide](./agents.md) | Consuming gapit from autonomous agents: schemas, introspection, errors |
| [FAQ](./faq.md) | Common questions, abricate differences, troubleshooting |

## Commands

| Command | What it does | Docs |
|---|---|---|
| `gapit screen` | Screen contig files or FASTQ reads for known genes | [Screening](./screen.md) |
| `gapit typing` | Designate phenotypes from screen result table(s) (typed gene databases) | [Typing schemes](./typing.md) |
| `gapit summary` | Summarize report table(s) into a gene presence/absence matrix | [Summary](./summary.md) |
| `gapit db fetch` | Fetch and build database(s) into the datadir (kaptive cluster dbs download on fetch) | [Databases](./databases.md) |
| `gapit db list` | List known databases (NAME, upstream PROVIDER, install state) | [Databases](./databases.md) |
| `gapit db install` | Install a local file after verifying its SHA256 | [Databases](./databases.md) |
| `gapit db build` | Build a database from a gene FASTA, or from GBK/GFF loci as a cluster database (optionally with a `--typing` phenotype spec) | [Custom databases](./custom-db.md) |
| `gapit setupdb` | Build BLAST indices for all databases under the datadir | here |
| `gapit schema` | Print the JSON Schema of a gapit output document (`report`, `typing_result`, `reads`, `reads2`, `cluster`, `summary`, `error`, `version`, `features`, `typing`) | [Outputs](./outputs.md) |
| `gapit mcp` | Run the MCP stdio server (also installed as the `gapit-mcp` console script) | [MCP server](./mcp.md) |

## Project documents

| File | Contents |
|---|---|
| [README](https://github.com/indexofire/gapit/blob/main/README.md) | Project overview, install, quick start, output contract, database table |
| [SPEC.md](https://github.com/indexofire/gapit/blob/main/SPEC.md) | The abricate behavior spec; source of truth for parity |
| [PLAN.md](https://github.com/indexofire/gapit/blob/main/PLAN.md) | Development roadmap, phase by phase |
| [AGENTS.md](https://github.com/indexofire/gapit/blob/main/AGENTS.md) | Contributor and agent guide: layout, conventions, verification |
| [CHANGELOG.md](https://github.com/indexofire/gapit/blob/main/CHANGELOG.md) | Change history |
