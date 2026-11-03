# gapit

[English](README.md) | [中文](README.zh.md)

Mass screening of contigs and reads for known genes — AMR, virulence, serotype, plasmid
replicons, species markers, or any custom database. Contig screening emits
abricate-compatible TSV; JSON and Markdown are first-class outputs. Python 3.12.

## Install

```bash
pip install gapit
```

BLAST+ (`blastn`, `blastx`, `makeblastdb`, `blastdbcmd`) and `minimap2` must be on
`PATH`. The pixi/conda paths install them for you; see
[installation](docs/installation.md) for those, Docker, and pinned setups.

```bash
# from a clone, with the external binaries managed for you
git clone https://github.com/indexofire/gapit.git
cd gapit
pixi install
pixi run gapit --version
```

## Quick start

Bundled databases work offline on first use; every other database is fetched with
`gapit db fetch`.

```bash
# contigs -> abricate-compatible TSV on stdout
gapit screen contigs.fa --db eco_dec

# JSON / Markdown for agents and humans
gapit screen contigs.fa --db eco_dec --format json
gapit screen contigs.fa --db eco_dec --format md

# FASTQ reads through minimap2
gapit screen --r1 R1.fastq.gz --r2 R2.fastq.gz --read-type sr --db eco_dec

# two-stage gene typing: screen writes the table, typing designates from it
gapit screen isolate.fa --db eco_dec --output result.tsv
gapit typing result.tsv

# summarize many reports into a gene matrix
gapit summary *.tsv

# databases
gapit db fetch -d ncbi
gapit db list

# introspect the output contract
gapit schema report
```

Contig defaults are `--minid 80 --mincov 80`; reads presence defaults to 90% alignment
breadth. Full flags and examples: [contigs](docs/screen.md), [reads](docs/reads.md),
[typing](docs/typing.md), [summary](docs/summary.md).

## Documentation

| Topic | Page |
|---|---|
| Install, quickstart | [installation](docs/installation.md), [quickstart](docs/quickstart.md) |
| Contig screening, cluster databases | [screen](docs/screen.md) |
| Reads and assembly surveys | [reads](docs/reads.md) |
| Gene matrices | [summary](docs/summary.md) |
| Typing schemes and phenotypes | [typing](docs/typing.md) |
| JSON/Markdown contracts, error envelope | [outputs](docs/outputs.md), `gapit schema` |
| Database catalog, providers, bundles | [databases](docs/databases.md) |
| Datadir, sharing, Docker | [datadir](docs/datadir.md) |
| Custom databases | [custom-db](docs/custom-db.md) |
| Agents and MCP server | [agents](docs/agents.md), [mcp](docs/mcp.md) |
| Optional gseek engine | [gseek](docs/gseek-engine.md) |

Full manual: <https://indexofire.github.io/gapit/>.

## Development and license

Contributor workflow and gates: [CONTRIBUTING.md](CONTRIBUTING.md); internal rules:
[AGENTS.md](AGENTS.md); trust boundaries: [SECURITY.md](SECURITY.md); changes:
[CHANGELOG.md](CHANGELOG.md).

gapit is MIT-licensed. It is a behavioral reimplementation of abricate (GPL-2.0) and copies
no Perl code. The seven gapit-curated bundled databases are permissively licensed; every
other provider downloads from upstream at fetch time under its own license. If you use
gapit, cite [CITATION.cff](CITATION.cff).
