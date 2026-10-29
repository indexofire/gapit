# gapit

Mass screening of contigs and reads for known genes — AMR, virulence, serotype, plasmid
replicons, species targets, or any custom database. An agent-first Python reimplementation
of [abricate](https://github.com/tseemann/abricate): byte-compatible TSV plus first-class
JSON and Markdown.

## Why gapit

- **Frozen abricate-compatible baseline.** Contig screening keeps the same BLAST
  pipeline, hit rules, and abricate-format TSV as a frozen baseline surface for existing
  pipelines. From v0.5.0 gapit pursues its own contract; abricate parity is no longer a
  release gate — the parity harness (`pixi run -e parity parity`) remains a regression
  reference for that frozen surface.
- **Machine-readable by design.** Versioned output schemas (`gapit.report/1`,
  `gapit.reads/1`, `gapit.summary/1`), self-describing via `gapit schema`, typed JSON error
  envelopes on stderr, documented exit codes. An agent can discover the whole contract
  without reading docs.
- **Reads, not just contigs.** `gapit screen --r1/--r2` screens FASTQ through minimap2, and
  `--r1` accepts assembly FASTA directly (content-detected, `map-ont` forced) for a fast
  presence survey. abricate cannot screen raw reads.

## Install

From PyPI (external BLAST+ and minimap2 binaries required on PATH — see
[docs/installation.md](docs/installation.md)):

```bash
pip install gapit
```

Or from a git clone with [pixi](https://pixi.sh), which manages the environment including
the external binaries:

```bash
git clone https://github.com/indexofire/gapit.git
cd gapit
pixi install
pixi run gapit --version
```

Alternatively, the repository's `environment.yml` puts gapit and the external binaries in a
conda environment:

```bash
git clone https://github.com/indexofire/gapit.git
cd gapit
conda env create -f environment.yml
conda activate gapit-env
gapit --version
```

The package supports Python 3.12 (3.12 only); the dev environment pins 3.12.

## Quick start

```bash
# Screen contigs (abricate-compatible TSV on stdout; --db is required)
gapit screen contigs.fa --db ncbi

# Agent- and human-readable outputs
gapit screen contigs.fa --db ncbi --format json
gapit screen contigs.fa --db ncbi --format md

# Two-stage designation: screen writes the gene table, typing designates from it
# (typed gene databases — the bundled ecoli_dec carries the dual-scheme DEC panel)
gapit screen isolates.fa --db ecoli_dec --output dec.tsv
gapit typing dec.tsv

# Screen FASTQ reads; --read-type picks the minimap2 preset
gapit screen --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --read-type sr --db card
gapit screen --r1 ont_reads.fastq.gz --db ncbi --read-type map-ont --format json

# Summarize report tables into a gene presence/absence matrix
gapit summary *.tsv

# Databases: eight ship in-wheel and work instantly; the rest download from upstream on fetch
gapit db fetch all         # default set (card, vfdb)
gapit db fetch ncbi
gapit db list
gapit db outdated                # flag databases older than the staleness threshold
gapit db search "tet(M)"         # look up genes across every installed database
gapit db build mydb my_genes.fa --tsv my_meta.tsv   # custom db from any FASTA

# Introspection
gapit schema report     # JSON Schema of gapit.report/1
```

Contig screening follows abricate defaults (`--minid 80`, `--mincov 80`). Reads-mode
presence defaults to 90% alignment breadth (`--min-breadth 90`).

## Databases

A database is a directory under the datadir, resolved from `--datadir`, then
`$GAPIT_DATADIR`, then `$XDG_DATA_HOME/gapit/db`, then `~/.local/share/gapit/db`
(use `gapit db path` to inspect it). On a fresh install
`gapit db list` bootstraps the datadir and lists the catalog (bundled rows ready for
first use); `gapit setupdb` creates it and materializes all eight bundled databases.
Twenty-three databases are
built into the catalog — twelve gene providers, seven kaptive cluster providers, and
eight bundled snapshots (`cje_cps`, `ecoh`, `ecoli_dec`, `lm_doumith`, `ncbi`, `resfinder`,
`vc_serogroup`, `vp_ok`; `ecoh`/`ncbi`/`resfinder` are also provider names). Every
provider builds on `gapit db fetch`,
which downloads from upstream — and eight content-provenance-audited ones also ship inside
the package (`ecoli_dec`, `lm_doumith`, and `vc_serogroup` plus the `ncbi`, `resfinder`,
`ecoh` snapshots — public domain / Apache-2.0 /
BSD-3-Clause / MIT — and the MIT typed-cluster snapshot `vp_ok`), materialized into the
datadir on first use with zero network; the
bundled `ecoli_dec` additionally carries per-gene identity floors (`floors.json`, `pic`
at 90%) that both engines apply before presence calls, and `cje_cps` ships floors plus
per-database `screen.json` defaults (MAPQ 20, depth 5, 10% lpxA relative depth). The
same `screen.json` (`gapit.screen/1`) mechanism lets any database declare its own
reads/contigs thresholds, so one global parameter set is not forced onto panels of
different amplicon lengths. Several
upstream licenses (CARD's McMaster terms, VFDB's CC BY-NC, Kaptive's GPL-3.0) forbid
redistribution inside an MIT-licensed distribution, so everything else is fetch-on-demand
— and `db fetch <name>` refreshes the provider-backed bundled names (`ecoh`, `ncbi`,
`resfinder`) from the latest upstream. For the database anatomy, datadir sharing,
conda/Docker guidance, and security notes, see [Data directory](docs/datadir.md).

| Name | Content | dbtype |
|---|---|---|
| `ncbi` | NCBI AMRFinderPlus curated AMR | nucl |
| `card` | CARD protein homolog resistance models | nucl |
| `resfinder` | CGE ResFinder acquired resistance genes | nucl |
| `argannot` | ARG-ANNOT acquired resistance genes | nucl |
| `plasmidfinder` | CGE PlasmidFinder replicons | nucl |
| `megares` | MEGARes antimicrobial resistance genes | nucl |
| `ecoh` | E. coli O and H antigens (srst2 EcOH) | nucl |
| `ecoli_dec` | Diarrheagenic E. coli marker panel (GB 4789.6 + risk-monitoring designation) | nucl |
| `vfdb` | VFDB virulence factors (set A, nucleotide) | nucl |
| `ecoli_vf` | E. coli virulence factors | nucl |
| `lm_doumith` | Listeria monocytogenes serogrouping (Doumith 2004) | nucl |
| `bacmet2` | BacMet2 biocide/resistance genes (protein) | prot |
| `victors` | Victors virulence factors | nucl |
| `upec_expec_vf` | UPEC/ExPEC virulence genes | nucl |
| `vc_serogroup` | V. cholerae O1/O139 serogroup + ctxA toxigenicity markers | nucl |
| `vp_ok` | V. parahaemolyticus O and K antigen loci (OAgc/CPSgc, vpasero; cluster, MIT, one-scan O:K serotype) | nucl |

gapit also reads datadirs built by abricate itself (legacy `~~~` headers). The reverse
does not hold: gapit-native databases use the `gapit/v1` header format (see SPEC.md §11),
which abricate cannot read. Protein databases such as `bacmet2` screen through blastx.

Gene **cluster** databases are a second kind: `gapit db build NAME loci.gbk|gff3` builds a
`kind: cluster` db (locus calls via minimap2, optional `--typing FILE` phenotype scoring),
and the seven Kaptive cluster databases (kpsc_k, kpsc_o, kosc_k, kosc_o, ab_k, ab_o,
ecoli_kps) fetch the Kaptive antigen-locus references at install time (GPL-3.0 content, so
nothing is bundled — cite Kaptive/Wyres et al. 2020 for results). The typed cluster
bundle `vp_ok`
merges the vpasero OAgc + CPSgc references (MIT) into one bundled db — a pure built-in in
the ecoli_dec/lm_doumith pattern, no `db fetch` — where a single
`gapit screen --db vp_ok` renders gapit.cluster/2 with `o_group`, `k_group`, and the
composed `serotype` O:K call (RIMD 2210633 → O3:K6), working straight after
`pip install gapit` with zero network.
`gapit db list --json`
surfaces each provider's content license where one is pinned.

## Output contract

- **stdout purity.** Data on stdout, diagnostics on stderr, always. `--quiet` silences
  stderr only.
- **Deterministic.** Stable sort orders, fixed tool parameters, no wall-clock timestamps
  inside data payloads.
- **Errors** print a JSON envelope to stderr and exit nonzero:

```text
{"schema": "gapit.error/1", "code": "...", "message": "...", "context": {...}}
```

| Exit code | Meaning |
|---|---|
| 0 | success |
| 1 | unexpected error |
| 2 | usage error |
| 3 | missing dependency |
| 4 | database error |
| 5 | input error |

## MCP server

gapit ships an MCP (Model Context Protocol) stdio server so agent runtimes can
screen assemblies without parsing CLI output: `gapit mcp` or the `gapit-mcp`
console script speaks newline-delimited JSON-RPC 2.0 on stdin/stdout (no extra
dependencies — the protocol is hand-rolled). It exposes ten tools:
`screen` (gapit.report/1 by default; `aligner minimap2` for a fast assembly
survey), `screen_reads` (FASTQ via minimap2, gapit.reads/1), `summary`
(gapit.summary/1), `schema`, `db_list`, plus the database tools `db_fetch`,
`db_build`, `db_search`, and `db_outdated` so an agent can self-provision
and inspect databases mid-session.
Tool failures return `isError: true` with the `gapit.error/1` envelope as text.
Register it with an MCP client:

```json
{"mcpServers": {"gapit": {"command": "gapit-mcp"}}}
```

## Optional gseek engine

`GAPIT_ENGINE=legacy|gseek|shadow` can route contig BLAST, protein `blastx-fast`,
reads minimap2, and cluster minimap2 through
[gseek](https://github.com/indexofire/gseek). `legacy` remains the default and
is retained on every path; `shadow` runs both and returns legacy. gseek requires
Python 3.12 and is installed by the pixi `default` environment (the `py312` feature) or
`pip install "gapit[gseek]"`. See [gseek engine](docs/gseek-engine.md).

## Development

| Task | Runs |
|---|---|
| `pixi run lint` | ruff check |
| `pixi run fmt` | ruff format |
| `pixi run fmt-check` | ruff format --check |
| `pixi run typecheck` | basedpyright (strict) |
| `pixi run test` | pytest, offline unit/CLI suite |
| `pixi run -e parity parity` | byte-diff screening vs real abricate |
| `pixi run -e parity summary-parity` | byte-diff summary vs real abricate |

User documentation lives in [`docs/index.md`](docs/index.md), rendered at
<https://indexofire.github.io/gapit/>.

`SPEC.md` is the parity contract, `AGENTS.md`/`CONTRIBUTING.md` the contributor rules,
`SECURITY.md` the trust-boundary and reporting policy, `CITATION.cff` the citation
metadata, and `CHANGELOG.md` the change history. `PLAN.md` is kept as a historical
roadmap (the current architecture lives in `AGENTS.md`, `SPEC.md`, and `docs/datadir.md`).

## License

gapit is MIT-licensed. It is a behavioral reimplementation of abricate (GPL-2.0) and
copies no Perl code; abricate itself remains GPL-2.0. The only database content inside
the package is the eight audited, permissively licensed bundles under `src/gapit/data/dbs/`
(public domain / Apache-2.0 / BSD-3-Clause / MIT); every other provider downloads from
upstream at fetch time, under its own license (SPEC.md §9).
