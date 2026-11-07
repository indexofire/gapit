# AGENTS.md — gspec

> Python tool for mass screening of contigs and reads for known genes, with gene-cluster
> typing and phenotype designation. An independent reimplementation of
> [abricate](https://github.com/tseemann/abricate); **from v0.5.0 onward gspec pursues its
> own contract** — abricate parity is a regression reference for the legacy table, not a
> design goal. **Agent-first**: outputs are machine-readable (JSON / Markdown) by design.

## 1. Mission

`gspec` answers one question: **which known genes are present in this assembly, and how confident
are we?** A modern, typed, testable Python tool that:

1. Ships the abricate-compatible 15-column contig TSV as a frozen baseline surface for existing
   pipelines (`SPEC.md` §5), while new development follows gspec's own design.
2. Adds first-class **JSON** and **Markdown** outputs so LLM agents and humans can consume
   results without parsing tab-delimited text.
3. Exposes stable, versioned output contracts (schemas, exit codes, error envelopes) that
   autonomous agents can rely on.

## 2. Environment & toolchain

- **Environment manager**: [pixi](https://pixi.sh) (conda-forge channel). Never use pip/conda
  directly; add dependencies to `pixi.toml` (`pixi add <pkg>` for conda, `pixi add --pypi <pkg>`
  for PyPI).
- **Python**: 3.12 in the pixi dev env; the package declares
  `requires-python = ">=3.12,<3.13"` and CI tests Python 3.12.
- **External binaries**: BLAST+ (`blastn`, `blastx`, `makeblastdb`, `blastdbcmd`) and `minimap2`,
  all from conda-forge/bioconda. Invoked only via `subprocess` with an argument list — never
  `shell=True`. Input normalization (fa/fq/gbk/embl, gz/bz2) is native (`seqconvert.py`);
  abricate's own `any2fasta` exists only inside the opt-in `parity` env, for the oracle.
- **Core libraries**: `click` (CLI), `pydantic` v2 (data models / JSON schema), `rich` (terminal
  output). No biopython — FASTA I/O is a small streaming parser we own.
- **Optional engine**: `gseek` (`gspec[gseek]`, Python 3.12) behind `GSPEC_ENGINE`;
  `legacy` remains the default and needs no gseek import.
- **Quality gates**: `ruff` (lint + format), `basedpyright` (strict mode), `pytest`.

### Commands (pixi tasks)

```bash
pixi run lint        # ruff check
pixi run fmt         # ruff format
pixi run fmt-check   # ruff format --check
pixi run typecheck   # basedpyright --strict
pixi run test        # pytest (unit, offline)
pixi run gspec       # the CLI itself
pixi run -e docs mkdocs build --strict     # bilingual docs site
pixi run -e parity parity                  # byte-diff screening vs real abricate
pixi run -e parity summary-parity          # byte-diff summary vs real abricate
```

Every change must leave `lint`, `fmt-check`, `typecheck`, and `test` green; doc changes also need
the strict docs build. Parity is an opt-in regression check for the legacy surface only.

## 3. Repository layout

A grouped map, not an inventory — open a directory for its full file list.

```
gspec/
├── AGENTS.md, CONTRIBUTING.md, SECURITY.md, CITATION.cff
├── SPEC.md              # frozen legacy table/header interop contract
├── PLAN.md              # historical development roadmap
├── docs/                # user docs (bilingual) + mkdocs nav
├── pixi.toml, pyproject.toml, recipe/meta.yaml
├── .github/workflows/ci.yml   # gates + scheduled clean-gates + parity job
├── src/gspec/
│   ├── cli.py, cmd_*.py       # click entrypoints: screen / summary / typing / db / schema / mcp
│   ├── config.py, dispatch.py, errors.py, proctools.py
│   ├── db.py, dbcodec.py, db_txn.py, records.py, dbbuild.py
│   ├── db_*_ops.py, bundled*.py, providers/   # fetch/build/rename/info/query + bundles
│   ├── blast.py, hits.py, screening*.py, fragments.py, minimap*.py, reads.py, readpairs.py
│   ├── cluster*.py, gbfeatures.py, gffparse.py, locus_map*.py
│   ├── gene_families.py, serotype_*.py, typing_*.py, gene_floors.py
│   ├── formats/               # tsv/csv/json/md/summary/cluster/typing_result/schemas/locus_map
│   ├── mcp*.py, gseek_engine.py
│   └── data/dbs/              # curated bundled databases (see §6)
├── scripts/                   # reproducible asset builders (shigella/, vibrio/) + calibration
└── tests/                     # data/, golden/, parity/, test_*.py
```

One file, one responsibility. New modules stay small (~250 LOC); the established CLI and
orchestration modules (`cmd_screen.py`, `db_ops.py`, `gseek_engine.py`) are deliberately larger.

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
- **Bilingual docs.** A user-facing doc change updates both `docs/foo.md` and `docs/foo.zh.md`.
- Match the style of the file you are editing; when in doubt, `ruff format` decides.

## 5. The agent-facing output contract (design center)

Treat this as a public API. Formats are detailed in `docs/outputs.md`, agent usage in
`docs/agents.md`, and every JSON document is introspectable with `gspec schema`.

- **Formats**: `--format tsv|csv|json|md` (default `tsv` for existing pipelines).
- **JSON**: every document carries a top-level `"schema": "<name>/<version>"`; `gspec schema
  --help` lists every registered document. Keys are snake_case and units explicit
  (`identity_pct`, `coverage_pct`). Semver the schema; never rename or retype a field in a
  minor bump.
- **Markdown**: YAML frontmatter (tool version, db, params, ISO-8601 UTC timestamp) + tables a
  human can read and an agent can regex reliably.
- **Errors**: failures print a JSON envelope to stderr
  `{"schema": "gspec.error/1", "code": "...", "message": "...", "context": {...}}` and exit with a
  documented code (1 unexpected, 2 usage, 3 missing dependency, 4 db error, 5 input error).
- **Typing** `[gspec-extension]`: two-stage designation — `gspec screen -o result.tsv --db NAME`
  detects genes (screen output is phenotype-free; typed and untyped gene dbs screen
  byte-identically), then `gspec typing result.tsv [-D|-f|-o|-q]` designates from the screen
  table (`gspec.typing_result/1`; cluster dbs refuse — their typing rides the screen).
  The MCP `typing` tool accepts the screen table text and returns the same document.
- **Databases** `[gspec-extension]`: `gspec db fetch|path|list|info|search|outdated|build|install|rename|setup`.
  Providers download from upstream at fetch time; only the seven gspec-curated bundles under
  `src/gspec/data/dbs/` ship inside the wheel (policy in §6). `db build NAME` with no input
  re-projects an installed gene database from its edited `records.jsonl` (the curation path);
  `db rename OLD NEW` is the migration path for names that predate the convention. Anatomy,
  sharing, and Docker policy live in `docs/databases.md` and `docs/datadir.md`.
  **Naming**: species databases use the genus initial plus two species letters (`cje`, `eco`,
  `vch`); genus-level databases use the first three genus letters (`shi`, `cam`), followed by a
  content suffix (`cje_cps`, `shi_serotype`). Published names stay stable.
  **Datadir precedence**: `--datadir`, `$GSPEC_DATADIR`, `$XDG_DATA_HOME/gspec/db` (absolute XDG
  only), then `~/.local/share/gspec/db`; `gspec db path [--json]` introspects it without
  creating, while `db list`/`db setup` bootstrap a missing root.
- **MCP** `[gspec-extension]`: `gspec mcp` / `gspec-mcp` stdio server exposing twelve tools:
  read-only `screen` (incl. `aligner minimap2`), `screen_reads`, `summary`, `typing`, `schema`,
  `db_path`, `db_search`, `db_outdated`; `db_list` may create the empty datadir on first use;
  the datadir-mutating `db_setup`, `db_fetch`, and `db_build` complete the set. Failures carry
  the `gspec.error/1` envelope.
- **stdout purity**: data on stdout, diagnostics on stderr, always. `--quiet` only affects stderr.
- **Self-description**: `gspec --version --json`, `gspec db list --json`, `gspec schema --help` —
  an agent must be able to discover everything without reading docs.

## 6. Bundled-database policy (what may ride the wheel)

`src/gspec/data/dbs/` is the sanctioned home for license-clean content shipped inside the
wheel: a directory qualifies by carrying `sequences` + `bundled.json` (auto-discovered by
`gspec.bundled`; `kind: cluster` in the metadata additionally requires `features.json`).
Materialization is kind-aware: gene bundles replay the gene-build pipeline, cluster
bundles the cluster-build pipeline (`bundled_cluster.py` rebuilds from the committed
locus FASTA + resolved feature table — gene ids stay snapshot-resolved, the typing spec
re-validates against them).
Bundling stays the exception, not the pattern: every provider fetches from upstream at
`db fetch` time. Third-party provider snapshots are never bundled, even under a
permissive license — the wheel carries gspec-curated content only. Every new bundle
passes three gates, in order:

1. **Curation gate** (hard): the content must be gspec-owned or gspec-curated from
   public-domain sources. A snapshot of a third-party provider database (ncbi, card,
   ecoh, resfinder, ...) does not qualify, regardless of its license.
2. **License gate** (hard, independent of location): only public-domain or permissive
   content (PD, Apache-2.0, BSD, MIT), provenance-audited. GPL (Kaptive) and
   non-commercial terms (CARD, VFDB) never ride the MIT distribution — moving files
   cannot fix this.
3. **Size gate** (budget): the bundled-root total stays ≤ 30 MiB, mechanically enforced by
   `tests/test_bundled.py` (every `pip install gspec` re-downloads all bundles; PyPI caps
   one file at 100 MiB).

Growth tiers past the budget — split rather than inflate the ceiling:

- **≤ 30 MiB**: keep adding to `src/gspec/data/dbs/`. Atomic data-with-code versioning and
  zero packaging machinery.
- **30–100 MiB**: split a separate `gspec-data` wheel with its own version (the
  `spacy[lookups]` → `spacy-lookups-data` pattern), pulled via a `gspec[datasets]` extra;
  `bundled_databases()` grows a second discovery root.
- **> 100 MiB, or fast-moving upstreams**: fetch-on-demand only. `db fetch` is the
  freshness path for every bundled name anyway — snapshots are point-in-time, and each
  bundled upstream refresh would otherwise force a wheel release.

## 7. Testing strategy

- **Unit**: pure functions (coverage %, merge rules, header parsing) — no I/O beyond `tests/data`.
- **Golden files**: fixed tiny db + fixed contigs → expected TSV/JSON/MD committed; update the
  committed files deliberately and review diffs like code.
- **Parity harness**: `pixi run -e parity parity` (and `summary-parity`) runs real abricate
  (conda) and gspec over a small genome corpus and diffs the gene calls byte-for-byte (file,
  gene, %identity, %coverage). **It is a regression reference for the frozen legacy surface
  (`SPEC.md` §2–§6), not a release gate for new work** — new features need not check abricate
  and may diverge deliberately.
- **Docs**: `pixi run -e docs mkdocs build --strict` must be warning-free; `docs/*.md` and
  `docs/*.zh.md` stay in sync.
- Full suite + CI on Python 3.12; parity is opt-in via the `parity` pixi env. Offline tests stay
  fast.

## 8. Domain knowledge (frozen arithmetic in SPEC.md; current contracts in docs/)

- **Databases** are directories `<datadir>/<dbname>/` holding `sequences` plus its BLAST index.
  Legacy headers are `>DB~~~GENE~~~ACCESSION~~~RESISTANCE` followed by a space and the product
  (`SPEC.md` §2); gspec-native headers use the `gspec/v1` codec whose `func` key carries a
  locked per-provider vocabulary (`docs/databases.md`). abricate cannot read gspec-native
  databases.
- **Legacy screening** = native normalize (`seqconvert.py`) → `blastn` against one db → 15-field
  rows → identity/coverage filter → **dedup identical `(contig, qstart, qend)`** (first row
  wins) → one TSV row per hit. The default path does not merge overlapping intervals;
  `--merge-fragments` is the explicit opt-in extension in `fragments.py`. Field values and the
  minimap arithmetic are frozen in `SPEC.md` §3–§5 — do not change them silently.
- **Cluster databases** are GBK/GFF-built locus panels screened by `minimap2 -x asm20` on
  assembly FASTA only. `--r1/--r2` and `--aligner minimap2` are usage errors against a cluster
  database: locus thresholds are assembly-calibrated, and read mapping to repeat-rich
  multi-kilobase loci needs a read preset, depth/MAPQ filtering, and truth-set calibration.
  Cluster builds (and `db setup`) write `sequences.asm20.mmi`; see `docs/screen.md`.
- **Reads** use streaming `minimap2` PAF with breadth/depth presence rules — `docs/reads.md`.
- **Typing** is a second stage for gene databases and rides the screen for cluster databases;
  rules, schemes, and output live in `docs/typing.md` and `gspec.typing/1|2`.

## 9. Git & workflow

- Conventional Commits (`feat:`, `fix:`, `test:`, `docs:`, `refactor:`).
- Never commit unless the user explicitly asks. Never force-push.
- `PLAN.md` is the historical roadmap; update it when scope changes, not retroactively.
  `CHANGELOG.md` records user-visible changes.

## 10. For agents working in this repo

1. Read `SPEC.md` before touching `blast.py`, `hits.py`, or `minimap.py`; read the linked doc
   before changing any other public surface.
2. Public contract changes (JSON schema, exit codes, CLI flags) require a schema version bump and
   a `CHANGELOG.md` entry (plus `PLAN.md` when scope changes).
3. Do not add dependencies without checking this file's §2 first — the list is intentionally short.
4. Verification is part of the task: `pixi run lint && pixi run fmt-check && pixi run typecheck
   && pixi run test` before declaring done; doc changes also run the strict docs build.
5. **Custom databases built with `gspec db build` are not done until provenance is written into
   the manifest** (`license`, `note`, `source_urls`, `upstream_version`): `note` must carry the
   literature citation + source software/version + accessions + build date; `license` must state
   the upstream license as declared (`unspecified` when the source states none); and whenever the
   content derives from a publication/third-party tool without explicit redistribution rights,
   `note` must append the objection clause ("if any author or rights holder objects to the use or
   redistribution of this database, please contact indexofire@gmail.com and the content will be
   removed or adjusted promptly upon verification"). Verify with `gspec db info -d NAME`.
   Also write the `info.json` sidecar with a one-line `description` (vendor/references/caveats/usage
   as applicable) so `db list` stays brief while `db info` carries the detail.
   Full rule: `docs/custom-db.md` § "Provenance and license metadata".
