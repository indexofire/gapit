# gspec documentation

gspec screens contig files and FASTQ reads for known genes — AMR, virulence, serotype,
plasmid replicons, species targets, or any custom database. It is an agent-first Python
reimplementation of [abricate](https://github.com/tseemann/abricate): abricate-compatible
TSV on stdout, plus first-class JSON and Markdown outputs with versioned schemas.

## Pages

| Page | Contents |
|---|---|
| [Installation](./installation.md) | PyPI wheel, conda environment file, and pixi paths, external binaries, verification, database bootstrap, shell completions |
| [Quickstart](./quickstart.md) | A complete first session on the bundled test fixture, offline |
| [Screening contigs](./screen.md) | `gspec screen` on FASTA/GBK/EMBL inputs: thresholds, filters, formats |
| [Screening reads](./reads.md) | FASTQ and assembly FASTA through minimap2: `--r1`/`--r2`, presets, breadth-based presence, two-stage survey |
| [Summary](./summary.md) | `gspec summary`: gene presence/absence matrix across report tables |
| [Databases](./databases.md) | The database catalog, `gspec db fetch/list/install`, datadirs, native `gapit/v1` format, cluster databases (GBK/GFF, typing.json, kaptive) |
| [Data directory](./datadir.md) | Datadir resolution (`db path`), database anatomy, sharing across envs, conda/Docker guidance, security/reproducibility |
| [Custom databases](./custom-db.md) | `gspec db build` walkthrough: any FASTA to a screenable database, worked examples |
| [Typing schemes](./typing.md) | The two-stage designation pipeline (`gspec screen` then `gspec typing`) and `gapit.typing/2`: rule/scheme primitives, the six-scheme cookbook (Doumith, ShigaTyper, meningotype, VP O/K, cholerae, DEC), allele probe tricks, calibration |
| [Meningococcal serogroup typing](./nme_sero.md) | `nme_sero` database: serogroup + ctrA + sodC depth control (`gspec screen` + `gspec typing`; no MLST/allele typing) |
| [Outputs](./outputs.md) | TSV/JSON/Markdown formats, schemas, error envelopes, exit codes |
| [MCP server](./mcp.md) | `gspec mcp`: analysis + database tools for agent runtimes over stdio JSON-RPC |
| [Agent guide](./agents.md) | Consuming gspec from autonomous agents: schemas, introspection, errors |
| [FAQ](./faq.md) | Common questions, abricate differences, troubleshooting |

## Commands

| Command | What it does | Docs |
|---|---|---|
| `gspec screen` | Screen contig files or FASTQ reads for known genes | [Screening](./screen.md) |
| `gspec typing` | Designate phenotypes from screen result table(s) (typed gene databases) | [Typing schemes](./typing.md) |
| `gspec summary` | Summarize report table(s) into a gene presence/absence matrix | [Summary](./summary.md) |
| `gspec db fetch` | Fetch and build database(s) into the datadir (kaptive cluster dbs download on fetch) | [Databases](./databases.md) |
| `gspec db list` | List known databases (NAME, upstream PROVIDER, install state); `-n` prints names only | [Databases](./databases.md) |
| `gspec db path` | Print the resolved datadir (no creation) | [Data directory](./datadir.md) |
| `gspec db install` | Install a local file after verifying its SHA256 | [Databases](./databases.md) |
| `gspec db build` | Build a database from a gene FASTA, or from GBK/GFF loci as a cluster database (optionally with a `--typing` phenotype spec) | [Custom databases](./custom-db.md) |
| `gspec db setup` | Bootstrap/index databases (`-a` all, `-d NAME` one; named providers fetch when missing) | [Data directory](./datadir.md) |
| `gspec schema` | Print the JSON Schema of a gspec output document (`report`, `typing_result`, `reads`, `reads2`, `cluster`, `summary`, `error`, `version`, `dblist`, `dbsetup`, `dboutdated`, `datadir`, `features`, `typing`, `floors`, `screen`) | [Outputs](./outputs.md) |
| `gspec mcp` | Run the MCP stdio server (also installed as the `gspec-mcp` console script) | [MCP server](./mcp.md) |

## Project documents

| File | Contents |
|---|---|
| [README](https://github.com/indexofire/gspec/blob/main/README.md) | Project overview, install, quick start, documentation map |
| [SPEC.md](https://github.com/indexofire/gspec/blob/main/SPEC.md) | Legacy abricate-table/header interop contract; parity is regression only |
| [PLAN.md](https://github.com/indexofire/gspec/blob/main/PLAN.md) | Development roadmap, phase by phase |
| [AGENTS.md](https://github.com/indexofire/gspec/blob/main/AGENTS.md) | Contributor and agent guide: layout, conventions, verification |
| [CHANGELOG.md](https://github.com/indexofire/gspec/blob/main/CHANGELOG.md) | Change history |
