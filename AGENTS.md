# AGENTS.md — gapit

> Python reimplementation of [abricate](https://github.com/tseemann/abricate): mass screening of
> contigs for antimicrobial resistance and virulence genes. **Agent-first**: every output is
> machine-readable (JSON / Markdown) by design, not as an afterthought.

## 1. Mission

`gapit` answers one question: **which known genes are present in this assembly, and how confident
are we?** It replaces abricate (Perl) with a modern, typed, testable Python tool that:

1. Produces byte-compatible TSV with abricate (drop-in replacement for existing pipelines).
2. Adds first-class **JSON** and **Markdown** outputs so LLM agents and humans can consume
   results without parsing tab-delimited text.
3. Exposes stable, versioned output contracts (schemas, exit codes, error envelopes) that
   autonomous agents can rely on.

## 2. Environment & toolchain

- **Environment manager**: [pixi](https://pixi.sh) (conda-forge channel). Never use pip/conda
  directly; add dependencies to `pixi.toml` (`pixi add <pkg>` for conda, `pixi add --pypi <pkg>`
  for PyPI).
- **Python**: 3.14 in the pixi dev env (current stable); the package declares
  `requires-python = ">=3.11"` and CI tests 3.11 / 3.13 / 3.14.
- **External binaries**: BLAST+ (`blastn`, `blastx`, `makeblastdb`, `blastdbcmd`) and `minimap2`
  (FASTQ read screening, SPEC.md §10), all from conda-forge/bioconda. Invoked only via `subprocess`
  with an argument list — never `shell=True`. Input normalization (fa/fq/gbk/embl, gz/bz2) is
  native (`seqconvert.py`); `any2fasta` is no longer a gapit dependency anywhere — the
  differential-validation oracle binary arrives transitively via abricate in the opt-in
  `parity` pixi env (prepend `.pixi/envs/parity/bin` to PATH for the differential test).
- **Core libraries**: `typer` (CLI), `pydantic` v2 (data models / JSON schema), `rich` (terminal
  output). No biopython — FASTA I/O is a small streaming parser we own.
- **Quality gates**: `ruff` (lint + format), `basedpyright` (strict mode), `pytest`.

### Commands (pixi tasks)

```bash
pixi run lint        # ruff check
pixi run fmt         # ruff format
pixi run typecheck   # basedpyright --strict
pixi run test        # pytest (unit, offline)
pixi run gapit       # the CLI itself
pixi run -e parity parity          # byte-diff screening vs real abricate (abricate-only env)
pixi run -e parity summary-parity  # byte-diff summary vs real abricate
```

Every change must leave `lint`, `typecheck`, and `test` green.

## 3. Repository layout

```
gapit/
├── AGENTS.md            # this file
├── PLAN.md              # development roadmap (phase-gated)
├── SPEC.md              # distilled abricate behavior spec (source of truth for parity)
├── pixi.toml
├── recipe/
│   └── meta.yaml        # conda recipe (submission deferred)
├── .github/workflows/
│   └── ci.yml           # gates matrix 3.11/3.13/3.14 + parity job
├── src/gapit/
│   ├── __init__.py
│   ├── cli.py           # typer entrypoint: screen / summary / db / list / setupdb / schema / mcp
│   ├── config.py        # datadir resolution, defaults, env vars
│   ├── dispatch.py      # shared CLI dispatch (error envelope → exit codes) + --datadir option
│   ├── proctools.py     # external-tool plumbing: argv subprocess runner + stderr notes
│   ├── fasta.py         # streaming FASTA reader + shared gz/bz2 text opener
│   ├── seqconvert.py    # native input normalization: fa/fq/gbk/embl (±gz/bz2) → FASTA
│   ├── db.py            # database discovery, header parsing, makeblastdb wrapper
│   ├── dbcodec.py       # gapit/v1 tagged-header codec (percent-encoded ids)
│   ├── records.py       # records.jsonl truth store + gapit.manifest/1 provenance
│   ├── dbbuild.py       # deterministic native-db build pipeline with self-check
│   ├── db_ops.py        # db use-cases: provider fetch + list (shared CLI + MCP; no typer)
│   ├── db_query_ops.py  # db use-cases: search + outdated over installed DBs (shared CLI + MCP)
│   ├── db_build_ops.py  # db use-case: custom FASTA+TSV → native db build (shared CLI + MCP)
│   ├── blast.py         # blastn invocation + tabular output parsing
│   ├── hits.py          # Hit model, identity/coverage computation, filtering, dedup
│   ├── minimap.py       # COVERAGE_MAP construction (exact abricate arithmetic)
│   ├── minimap2_run.py  # minimap2 invocation layer for reads mode (streaming PAF, --cs/NM tags)
│   ├── paf.py           # PAF row parsing + interval arithmetic (minimap2 output boundary)
│   ├── report.py        # Report model: the canonical in-memory result
│   ├── screening.py     # blastn screen use-case + shared engine helpers (OutputFormat, AlignerEnum)
│   ├── screening_reads.py # minimap2 use-cases: --r1/--r2 reads + --aligner minimap2 assemblies
│   ├── reads.py         # FASTQ mode: minimap2 PAF parsing, coverage breadth/depth, presence
│   ├── summary.py       # summary core: parse report tables into a gene matrix
│   ├── cmd_screen.py    # `gapit screen` CLI (registered from cli.py)
│   ├── cmd_summary.py   # `gapit summary` CLI (registered from cli.py)
│   ├── cmd_db.py        # `gapit db` command group (fetch | list; registers the subcommands)
│   ├── cmd_db_install.py # `gapit db install`: SHA256-verified local-file install
│   ├── cmd_db_build.py  # `gapit db build` CLI (custom FASTA → native db)
│   ├── cmd_db_search.py # `gapit db search` CLI (records.jsonl lookup)
│   ├── cmd_db_outdated.py # `gapit db outdated` CLI (staleness report)
│   ├── mcp.py           # MCP stdio server (hand-rolled JSON-RPC 2.0); backs gapit-mcp
│   ├── mcp_tools.py     # MCP tool implementations (call the shared use-cases)
│   ├── mcp_schemas.py   # MCP tools/list declarations (names, descriptions, inputSchemas)
│   ├── errors.py        # typed errors + JSON error envelope
│   ├── formats/
│   │   ├── tsv.py       # abricate-compatible TSV/CSV
│   │   ├── json.py      # versioned JSON (gapit.report/1 et al.)
│   │   ├── md.py        # Markdown (human + agent readable, YAML frontmatter)
│   │   ├── schemas.py   # registered output models behind `gapit schema`
│   │   └── summary.py   # summary matrix renderers (TSV/CSV/JSON/MD)
│   ├── providers/       # 12 DB providers + common.py helpers + snapshots.py loader
│   ├── data/snapshots/  # bundled card + vfdb snapshot archives (.tar.gz)
│   └── py.typed
└── tests/
    ├── data/            # tiny synthetic db + contigs + reads (fast, offline)
    ├── golden/          # expected outputs incl. abricate reference TSVs
    ├── parity/          # corpus + run_parity.py / run_summary_parity.py (opt-in env)
    └── test_*.py        # unit + CLI + provider + integration tests
```

One file, one responsibility. Target ≤ 250 LOC per module; split before it hurts.

## 4. Coding conventions

- **Strict typing everywhere.** basedpyright strict; no `Any` unless isolated and justified in a
  comment. No `# type: ignore`, no `cast` to silence real errors.
- **Parse, don't validate.** BLAST rows, FASTA records, and db headers become typed models at the
  boundary; the core logic never touches raw strings.
- **No silent failures.** Errors are typed (`errors.py`), carry context, and map to documented
  exit codes. Empty `except` blocks are forbidden.
- **Deterministic output.** Stable sort orders, no wall-clock timestamps inside data payloads
  (metadata block only), LF line endings, UTF-8.
- **TDD for core logic.** hit filtering, merging, and coverage-map math are written test-first.
- Match the style of the file you are editing; when in doubt, `ruff format` decides.

## 5. The agent-facing output contract (design center)

This is what distinguishes gapit from abricate. Treat it as a public API.

- **Formats**: `--format tsv|csv|json|md` (default `tsv` for abricate compatibility).
- **JSON**: top-level `"schema": "gapit.report/1"`; schema introspectable via
  `gapit schema report | reads | reads2 | summary | list | error | version`. Keys are snake_case,
  units explicit (`identity_pct`, `coverage_pct`). Semver the schema; never rename or
  retype a field in a minor bump.
- **Markdown**: YAML frontmatter (tool version, db, params, ISO-8601 UTC timestamp) + tables a
  human can read and an agent can regex reliably.
- **Errors**: failures print a JSON envelope to stderr
  `{"schema": "gapit.error/1", "code": "...", "message": "...", "context": {...}}` and exit with a
  documented non-zero code (2 = usage, 3 = missing dependency, 4 = db error, 5 = input error).
- **DB acquisition** `[gapit-extension]`: `gapit db fetch|list|search|outdated|build|install` —
  provider fetch (bundled card/vfdb snapshots install offline; `--from-source` forces upstream),
  provider listing, records.jsonl gene search, staleness report, custom FASTA→native-db build,
  and SHA256-verified local-file install.
- **MCP** `[gapit-extension]`: `gapit mcp` / `gapit-mcp` stdio server exposing nine tools:
  read-only `screen` (incl. `aligner minimap2` assembly survey), `screen_reads` (FASTQ via
  minimap2), `summary`, `schema`, `db_list`, `db_search`, `db_outdated`, plus the datadir-mutating
  `db_fetch` (installs provider databases; may download) and `db_build` (writes a custom db);
  tool failures carry the `gapit.error/1` envelope.
- **stdout purity**: data on stdout, diagnostics on stderr, always. `--quiet` only affects stderr.
- **Self-description**: `gapit --version --json`, `gapit list --json`, `gapit schema` — an agent
  must be able to discover everything without reading docs.

## 6. Testing strategy

- **Unit**: pure functions (coverage %, merge rules, header parsing) — no I/O beyond `tests/data`.
- **Golden files**: fixed tiny db + fixed contigs → expected TSV/JSON/MD committed; update the
  committed files deliberately and review diffs like code.
- **Parity harness**: `pixi run -e parity parity` (and `summary-parity`) runs real abricate
  (conda) and gapit over a small genome corpus and diffs the gene calls byte-for-byte (file,
  gene, %identity, %coverage). Parity on the corpus is the release gate for v1.0.
- Full suite + CI matrix 3.11/3.13/3.14; parity byte-diff is opt-in via the `parity` pixi env.
  Offline tests stay fast.

## 7. Domain knowledge (abridged; full detail in SPEC.md)

- A **database** is a directory `<datadir>/<dbname>/sequences`: a nucleotide FASTA whose headers
  encode `>~~~GENE~~~PRODUCT` (and accession / resistance metadata depending on the source db).
  `makeblastdb -dbtype nucl` builds the index. gapit also writes its own native `gapit/v1`
  tagged-header format (SPEC §11), whose `func` key carries a locked per-provider function
  vocabulary (antibiotic classes, `virulence`, `replicon`, ...); abricate cannot read
  gapit-native databases.
- Screening = native normalize (`seqconvert.py`, any2fasta-equivalent semantics) → `blastn` of
  query contigs against one db → 15-field
  tabular hits → filter by identity / coverage thresholds → **dedup hits sharing identical
  `(contig, qstart, qend)`** (first/best BLAST row wins) → one TSV row per surviving hit.
  abricate does **not** merge overlapping intervals — do not "improve" this on the default path.
- Key computed fields: `%COVERAGE = 100*(length-gaps)/slen` (filtered unrounded, displayed
  `%.2f`), `%IDENTITY` (BLAST pident, never post-filtered), `COVERAGE_MAP` (15-char minimap),
  `GAPS`. Exact formulas, the dedup rule, and the minimap arithmetic live in `SPEC.md` — when
  abricate and intuition disagree, **abricate wins** (parity is a feature).

## 8. Git & workflow

- Conventional Commits (`feat:`, `fix:`, `test:`, `docs:`, `refactor:`).
- Never commit unless the user explicitly asks. Never force-push.
- `PLAN.md` tracks phases; update it when scope changes, not retroactively.

## 9. For agents working in this repo

1. Read `SPEC.md` before touching `blast.py`, `hits.py`, or `minimap.py`.
2. Public contract changes (JSON schema, exit codes, CLI flags) require a schema version bump and
   a note in `PLAN.md`.
3. Do not add dependencies without checking this file's §2 first — the list is intentionally short.
4. Verification is part of the task: `pixi run lint && pixi run typecheck && pixi run test` before
   declaring done.
