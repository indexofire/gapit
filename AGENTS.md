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
- **External binaries**: BLAST+ (`blastn`, `blastx`, `makeblastdb`, `blastdbcmd`),
  `any2fasta` (input normalization: gbk/embl/gz/bz2), and `minimap2` (FASTQ read screening,
  SPEC.md §10), all from conda-forge/bioconda. Invoked only via `subprocess` with an argument
  list — never `shell=True`.
- **Core libraries**: `typer` (CLI), `pydantic` v2 (data models / JSON schema), `rich` (terminal
  output). No biopython — FASTA I/O is a small streaming parser we own.
- **Quality gates**: `ruff` (lint + format), `basedpyright` (strict mode), `pytest`.

### Commands (pixi tasks)

```bash
pixi run lint        # ruff check
pixi run fmt         # ruff format
pixi run typecheck   # basedpyright --strict
pixi run test        # pytest (unit, offline)
pixi run parity      # golden-file diff against real abricate (requires abricate in env)
pixi run gapit       # the CLI itself
```

Every change must leave `lint`, `typecheck`, and `test` green.

## 3. Repository layout (target)

```
gapit/
├── AGENTS.md            # this file
├── PLAN.md              # development roadmap (phase-gated)
├── SPEC.md              # distilled abricate behavior spec (source of truth for parity)
├── pixi.toml
├── src/gapit/
│   ├── __init__.py
│   ├── cli.py           # typer entrypoint: screen / summary / setupdb / list / schema
│   ├── config.py        # datadir resolution, defaults, env vars
│   ├── fasta.py         # streaming FASTA reader (plain + gz)
│   ├── db.py            # database discovery, header parsing, makeblastdb wrapper
│   ├── blast.py         # blastn invocation + tabular output parsing
│   ├── hits.py          # Hit model, identity/coverage computation, filtering, dedup
│   ├── minimap.py       # COVERAGE_MAP construction (exact abricate arithmetic)
│   ├── report.py        # Report model: the canonical in-memory result
│   ├── screening.py     # the screen use-case: validation, orchestration, format dispatch
│   ├── reads.py         # FASTQ mode: minimap2 PAF parsing, coverage breadth/depth, presence
│   ├── summary.py       # multi-file summary matrix
│   ├── formats/
│   │   ├── tsv.py       # abricate-compatible TSV/CSV
│   │   ├── json.py      # versioned JSON (schema: gapit.report/1)
│   │   └── md.py        # Markdown (human + agent readable, YAML frontmatter)
│   ├── errors.py        # typed errors + JSON error envelope
│   └── py.typed
└── tests/
    ├── data/            # tiny synthetic db + contigs (fast, offline)
    ├── golden/          # expected outputs incl. abricate reference TSVs
    └── ...
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
  `gapit schema report | summary | error`. Keys are snake_case, units explicit (`identity_pct`,
  `coverage_pct`). Semver the schema; never rename or retype a field in a minor bump.
- **Markdown**: YAML frontmatter (tool version, db, params, ISO-8601 UTC timestamp) + tables a
  human can read and an agent can regex reliably.
- **Errors**: failures print a JSON envelope to stderr
  `{"schema": "gapit.error/1", "code": "...", "message": "...", "context": {...}}` and exit with a
  documented non-zero code (2 = usage, 3 = missing dependency, 4 = db error, 5 = input error).
- **stdout purity**: data on stdout, diagnostics on stderr, always. `--quiet` only affects stderr.
- **Self-description**: `gapit --version --json`, `gapit list --json`, `gapit schema` — an agent
  must be able to discover everything without reading docs.

## 6. Testing strategy

- **Unit**: pure functions (coverage %, merge rules, header parsing) — no I/O beyond `tests/data`.
- **Golden files**: fixed tiny db + fixed contigs → expected TSV/JSON/MD committed; regenerate via
  `pixi run golden --update`, review diffs like code.
- **Parity harness**: `pixi run parity` runs real abricate (conda) and gapit over a small genome
  corpus and diffs the gene calls (file, gene, %identity, %coverage). Parity on the corpus is the
  release gate for v1.0.
- Tests must run offline and fast (< 30 s) except `parity`, which is opt-in.

## 7. Domain knowledge (abridged; full detail in SPEC.md)

- A **database** is a directory `<datadir>/<dbname>/sequences`: a nucleotide FASTA whose headers
  encode `>~~~GENE~~~PRODUCT` (and accession / resistance metadata depending on the source db).
  `makeblastdb -dbtype nucl` builds the index.
- Screening = `any2fasta` normalize → `blastn` of query contigs against one db → 15-field
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
