# gspec

[English](README.md) | [中文](README.zh.md)

**Declarative typing databases for bacterial genomes.** One directory describes a
database (records, provenance, licence, optional sidecars); one versioned scheme
document decides the phenotype; BLAST, minimap2 or k-mer containment supplies the
evidence. AMR, virulence, serotype, antigen-locus and custom panels all speak the
same language. Contig screening emits abricate-compatible TSV; JSON and Markdown are
first-class outputs. Python 3.12.

> Formerly `gapit` (≤0.8.1). `gspec` 0.9 renamed every on-disk contract too;
> pre-0.9 databases and archived documents are still read. The `gapit`
> compatibility package (`pip install gapit`) keeps `import gapit` and the
> `gapit` command. See [Migrating from gapit](docs/migration.md).

## Install

```bash
pip install gspec
```

BLAST+ (`blastn`, `blastx`, `makeblastdb`, `blastdbcmd`) and `minimap2` must be on
`PATH`. The pixi/conda paths install them for you; see
[installation](docs/installation.md) for those, Docker, and pinned setups.

```bash
# from a clone, with the external binaries managed for you
git clone https://github.com/indexofire/gspec.git
cd gspec
pixi install
pixi run gspec --version
```

## Quick start

Bundled databases work offline on first use; every other database is fetched with
`gspec db fetch`.

```bash
# contigs -> abricate-compatible TSV on stdout
gspec screen contigs.fa --db eco_dec

# JSON / Markdown for agents and humans
gspec screen contigs.fa --db eco_dec --format json
gspec screen contigs.fa --db eco_dec --format md

# FASTQ reads through minimap2
gspec screen --r1 R1.fastq.gz --r2 R2.fastq.gz --read-type sr --db eco_dec

# two-stage gene typing: screen writes the table, typing designates from it
gspec screen isolate.fa --db eco_dec --output result.tsv
gspec typing result.tsv

# summarize many reports into a gene matrix
gspec summary *.tsv

# databases
gspec db fetch -d ncbi
gspec db list

# introspect the output contract
gspec schema report
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
| JSON/Markdown contracts, error envelope | [outputs](docs/outputs.md), `gspec schema` |
| Database catalog, providers, bundles | [databases](docs/databases.md) |
| Datadir, sharing, Docker | [datadir](docs/datadir.md) |
| Custom databases | [custom-db](docs/custom-db.md) |
| Agents and MCP server | [agents](docs/agents.md), [mcp](docs/mcp.md) |
| Optional gseek engine | [gseek](docs/gseek-engine.md) |

Full manual: <https://indexofire.github.io/gspec/>.

## Development and license

Contributor workflow and gates: [CONTRIBUTING.md](CONTRIBUTING.md); internal rules:
[AGENTS.md](AGENTS.md); trust boundaries: [SECURITY.md](SECURITY.md); changes:
[CHANGELOG.md](CHANGELOG.md).

gspec is MIT-licensed. It is a behavioral reimplementation of abricate (GPL-2.0) and copies
no Perl code. The seven gspec-curated bundled databases are permissively licensed; every
other provider downloads from upstream at fetch time under its own license. If you use
gspec, cite [CITATION.cff](CITATION.cff).
