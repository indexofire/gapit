# gapit

Mass screening of contigs and reads for antimicrobial resistance and virulence genes. An
agent-first Python reimplementation of [abricate](https://github.com/tseemann/abricate):
byte-compatible TSV plus first-class JSON and Markdown.

## Why gapit

- **Drop-in abricate replacement.** Same BLAST pipeline, same hit rules, abricate-format
  TSV on stdout. Parity with abricate 1.4.0 is a release gate, checked by diffing gene
  calls against real abricate (`pixi run -e parity parity`).
- **Machine-readable by design.** Versioned output schemas (`gapit.report/1`,
  `gapit.reads/1`, `gapit.summary/1`), self-describing via `gapit schema`, typed JSON error
  envelopes on stderr, documented exit codes. An agent can discover the whole contract
  without reading docs.
- **Reads, not just contigs.** `gapit screen --r1/--r2` screens FASTQ through minimap2, and
  `--r1` accepts assembly FASTA directly (content-detected, `map-ont` forced) for a fast
  presence survey. abricate cannot screen raw reads.

## Install

[pixi](https://pixi.sh) manages the environment, including the external binaries
(BLAST+, minimap2):

```bash
git clone https://github.com/indexofire/gapit.git
cd gapit
pixi install
pixi run gapit --version
```

The package supports Python 3.11+; the dev environment pins 3.14.

## Quick start

```bash
# Screen contigs (abricate-compatible TSV on stdout, --db defaults to ncbi)
gapit screen contigs.fa --db ncbi

# Agent- and human-readable outputs
gapit screen contigs.fa --db ncbi --format json
gapit screen contigs.fa --format md

# Screen FASTQ reads; --read-type picks the minimap2 preset
gapit screen --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --read-type sr --db card
gapit screen --r1 ont_reads.fastq.gz --read-type map-ont --format json

# Summarize report tables into a gene presence/absence matrix
gapit summary *.tsv

# Databases: card + vfdb install offline from bundled snapshots, others fetch on demand
gapit db fetch
gapit db fetch ncbi
gapit db list
gapit db outdated                # flag stale databases and newer bundled snapshots
gapit db search "tet(M)"         # look up genes across every installed database
gapit db build mydb my_genes.fa --tsv my_meta.tsv   # custom db from any FASTA

# Introspection
gapit list              # installed databases (abricate --list compatible)
gapit schema report     # JSON Schema of gapit.report/1
```

Contig screening follows abricate defaults (`--minid 80`, `--mincov 80`). Reads-mode
presence defaults to 90% alignment breadth (`--min-breadth 90`).

## Databases

A database is a directory under the datadir, resolved from `$GAPIT_DATADIR`, then
`~/.local/share/gapit/db` (override per call with `--datadir`). Twelve providers exist;
`card` and `vfdb` ship inside the wheel and install with zero network, the rest download
from upstream when fetched.

| Name | Content | dbtype |
|---|---|---|
| `ncbi` | NCBI AMRFinderPlus curated AMR (default db) | nucl |
| `card` | CARD protein homolog resistance models | nucl |
| `resfinder` | CGE ResFinder acquired resistance genes | nucl |
| `argannot` | ARG-ANNOT acquired resistance genes | nucl |
| `plasmidfinder` | CGE PlasmidFinder replicons | nucl |
| `megares` | MEGARes antimicrobial resistance genes | nucl |
| `ecoh` | E. coli O and H antigens (srst2 EcOH) | nucl |
| `vfdb` | VFDB virulence factors (set A, nucleotide) | nucl |
| `ecoli_vf` | E. coli virulence factors | nucl |
| `bacmet2` | BacMet2 biocide/resistance genes (protein) | prot |
| `victors` | Victors virulence factors | nucl |
| `upec_expec_vf` | UPEC/ExPEC virulence genes | nucl |

gapit also reads datadirs built by abricate itself (legacy `~~~` headers). The reverse
does not hold: gapit-native databases use the `gapit/v1` header format (see SPEC.md §11),
which abricate cannot read. Protein databases such as `bacmet2` screen through blastx.

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
dependencies — the protocol is hand-rolled). It exposes eight tools:
`screen` (gapit.report/1 by default), `summary` (gapit.summary/1), `schema`,
`db_list`, plus the database tools `db_fetch`, `db_build`, `db_search`, and
`db_outdated` so an agent can self-provision and inspect databases mid-session.
Tool failures return `isError: true` with the `gapit.error/1` envelope as text.
Register it with an MCP client:

```json
{"mcpServers": {"gapit": {"command": "gapit-mcp"}}}
```

## Development

| Task | Runs |
|---|---|
| `pixi run lint` | ruff check |
| `pixi run fmt` | ruff format |
| `pixi run typecheck` | basedpyright (strict) |
| `pixi run test` | pytest, 374 offline tests |
| `pixi run -e parity parity` | byte-diff screening vs real abricate |
| `pixi run -e parity summary-parity` | byte-diff summary vs real abricate |

User documentation lives in [`docs/index.md`](docs/index.md), rendered at
<https://indexofire.github.io/gapit/>.

`SPEC.md` is the parity contract, `PLAN.md` the roadmap, `AGENTS.md` the contributor
guide, `CHANGELOG.md` the change history.

## License

GPL-2.0-compatible: gapit is a behavioral reimplementation of GPL-2.0 abricate and copies
no Perl code. Bundled database content retains its original upstream licenses
(SPEC.md §9).
