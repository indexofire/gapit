# AGENTS.md — gapit

> Python tool for mass screening of contigs and reads for known genes, with gene-cluster
> typing and phenotype designation. Originally a reimplementation of
> [abricate](https://github.com/tseemann/abricate); **from v0.5.0 onward gapit pursues its
> own contract** — abricate parity is no longer a design goal. **Agent-first**: every
> output is machine-readable (JSON / Markdown) by design, not as an afterthought.

## 1. Mission

`gapit` answers one question: **which known genes are present in this assembly, and how confident
are we?** A modern, typed, testable Python tool that:

1. Ships an abricate-compatible contig TSV as the historical baseline surface (frozen
   behavior, kept for existing pipelines), while new development follows gapit's own design.
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
- **External binaries**: BLAST+ (`blastn`, `blastx`, `makeblastdb`, `blastdbcmd`) and `minimap2`
  (FASTQ read screening, SPEC.md §10), all from conda-forge/bioconda. Invoked only via `subprocess`
  with an argument list — never `shell=True`. Input normalization (fa/fq/gbk/embl, gz/bz2) is
  native (`seqconvert.py`); `any2fasta` is no longer a gapit dependency anywhere — the
  differential-validation oracle binary arrives transitively via abricate in the opt-in
  `parity` pixi env (prepend `.pixi/envs/parity/bin` to PATH for the differential test).
- **Core libraries**: `click` (CLI), `pydantic` v2 (data models / JSON schema), `rich` (terminal
  output). No biopython — FASTA I/O is a small streaming parser we own.
- **Optional engine**: `gseek` (`gapit[gseek]`, Python 3.12) behind `GAPIT_ENGINE`;
  `legacy` remains the default and needs no gseek import.
- **Quality gates**: `ruff` (lint + format), `basedpyright` (strict mode), `pytest`.

### Commands (pixi tasks)

```bash
pixi run lint        # ruff check
pixi run fmt         # ruff format
pixi run fmt-check   # ruff format --check
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
├── CONTRIBUTING.md      # contributor workflow and gates
├── SECURITY.md          # trust boundaries + vulnerability reporting
├── CITATION.cff         # citation metadata
├── PLAN.md              # historical development roadmap
├── SPEC.md              # frozen contract for the baseline surface (distilled from abricate)
├── pixi.toml
├── recipe/
│   └── meta.yaml        # conda recipe (submission deferred)
├── .github/workflows/
│   └── ci.yml           # gates matrix Python 3.12 + parity job
├── src/gapit/
│   ├── __init__.py
│   ├── cli.py           # click entrypoint: screen / summary / db / schema / mcp
│   ├── config.py        # datadir resolution, first-run bootstrap, defaults, env vars
│   ├── dispatch.py      # shared CLI dispatch (error envelope → exit codes) + --datadir option + --output target
│   ├── proctools.py     # external-tool plumbing: argv subprocess runner + stderr notes
│   ├── fasta.py         # streaming FASTA reader + shared gz/bz2 text opener
│   ├── seqconvert.py    # native input normalization: fa/fq/gbk/embl (±gz/bz2) → FASTA
│   ├── db.py            # database discovery, header parsing, makeblastdb wrapper
│   ├── dbcodec.py       # gapit/v1 tagged-header codec (percent-encoded ids)
│   ├── db_txn.py        # per-db write lock + staging build + rename publish
│   ├── records.py       # records.jsonl truth store + gapit.manifest/1 provenance
│   ├── dbbuild.py       # deterministic native-db build pipeline with self-check
│   ├── clusterbuild.py  # cluster-db build pipeline: GBK/GFF → locus FASTA + features.json + typing copy
│   ├── cluster.py       # cluster-screening engine core: minimap2 asm20 → per-locus/per-gene calls
│   ├── cluster_math.py  # cs-walk → union coverage/identity math + verdicts (pure functions)
│   ├── db_ops.py        # db use-cases: fetch + list (shared CLI + MCP; no typer)
│   ├── db_setup_ops.py  # db use-case: bootstrap + materialize + index (gapit.dbsetup/1)
│   ├── db_rename_ops.py # db use-case: rename an installed database (gene rebuild / cluster re-stamp)
│   ├── db_info_ops.py   # db use-case: provenance/literature/caveats document (gapit.dbinfo/1)
│   ├── db_info_manual.py # `db info` man-style renderer (sections, wrapping, terminal width)
│   ├── db_info_models.py # leaf models shared by bundled.json and db info (references/guidance)
│   ├── db_query_ops.py  # db use-cases: search + outdated over installed DBs (shared CLI + MCP)
│   ├── db_build_ops.py  # db use-case: custom FASTA+TSV → gene db, GBK/GFF → cluster db (shared CLI + MCP)
│   ├── db_build_meta.py # `db build --tsv` metadata sidecar: parse + merge (shared by the build use-case)
│   ├── gbfeatures.py    # GenBank FEATURES/ORIGIN parser → LocusFeatures/GeneFeature + gapit.features/1
│   ├── gene_floors.py   # gapit.floors/1 per-gene identity floors: model, loader, build validation, reads gate
│   ├── gffparse.py      # GFF3 parser (embedded ##FASTA or sidecar) → the same locus/gene models
│   ├── typing_models.py # gapit.typing schema (rules + named schemes; /1 degrades to one default scheme)
│   ├── typing_rules.py  # gapit.typing rule models (weighted_genes/exact_set/cluster_match/learned_linear) + reference checks
│   ├── typing_decide.py # typing decision layer: ranking, control-gene gate, unique-group mixed override
│   ├── typing_results.py # typing evaluation-result models (PhenotypeDetail, SchemeCall)
│   ├── typing_engine.py # typing evaluation, cluster path: rule scoring → phenotype call + breakdown
│   ├── typing_gene.py   # typing evaluation, gene path (DEC markers): hits → per-scheme report phenotypes
│   ├── typing_cluster2.py # typing evaluation, multi-scheme cluster path: report → per-scheme phenotypes map
│   ├── bundled.py       # bundled databases: discovery + kind-forked materialize-into-datadir
│   ├── bundled_cluster.py # cluster-bundle materialization: committed artifacts → cluster build
│   ├── blast.py         # blastn invocation + tabular output parsing
│   ├── hits.py          # Hit model, identity/coverage computation, filtering, dedup
│   ├── fragments.py     # opt-in --merge-fragments cross-contig gene-fragment merging
│   ├── minimap.py       # COVERAGE_MAP construction (exact abricate arithmetic)
│   ├── minimap2_run.py  # minimap2 invocation layer for reads mode (streaming PAF, --cs/NM tags)
│   ├── paf.py           # PAF row parsing + interval arithmetic (minimap2 output boundary)
│   ├── report.py        # Report model: the canonical in-memory result
│   ├── engines.py       # shared OutputFormat / AlignerEnum (leaf module for the three engines)
│   ├── gseek_engine.py  # optional gseek bridge (legacy/gseek/shadow) for BLAST/minimap2
│   ├── progress.py      # stderr-bound rich batch progress bar
│   ├── screening.py      # blastn screen use-case + the shared run_screen dispatcher
│   ├── screening_cluster.py # cluster screen use-case: kind guards, typing wiring, rendering
│   ├── screening_reads.py # minimap2 use-cases: --r1/--r2 reads + --aligner minimap2 assemblies
│   ├── screening_reads_positional.py # positional FASTQ wildcard: pairing + per-sample screens
│   ├── typing_input.py   # `gapit typing` input: contig/reads tables → per-FILE folded GeneCalls
│   ├── typing_table.py   # `gapit typing` table parsing: contig 15-column + reads #SAMPLE/BREADTH% shapes
│   ├── typing_ops.py     # `gapit typing` use-case: resolve db, run the engine, render
│   ├── reads.py          # FASTQ mode: minimap2 PAF parsing, coverage breadth/depth, presence
│   ├── readpairs.py      # positional FASTQ filename → sample/lane inference
│   ├── summary.py       # summary core: parse report tables into a gene matrix
│   ├── cmd_screen.py    # `gapit screen` CLI (registered from cli.py)
│   ├── cmd_screen_positionals.py # screen positional guard: db-name-as-file usage error
│   ├── cmd_screen_reads_args.py # reads-mode CLI arg rules (comma split + flag rejects)
│   ├── cmd_summary.py   # `gapit summary` CLI (registered from cli.py)
│   ├── cmd_typing.py    # `gapit typing` CLI (registered from cli.py)
│   ├── cmd_db.py        # `gapit db` command group (fetch/path/list/build/install/search/outdated/setup)
│   ├── cmd_db_install.py # `gapit db install`: SHA256-verified local-file install
│   ├── cmd_db_rename.py # `gapit db rename` CLI (naming-convention migration)
│   ├── cmd_db_info.py   # `gapit db info -d NAME` CLI (manual page / TSV / JSON)
│   ├── cmd_db_build.py  # `gapit db build` CLI (custom FASTA → native db)
│   ├── cmd_db_search.py # `gapit db search` CLI (records.jsonl lookup)
│   ├── cmd_db_outdated.py # `gapit db outdated` CLI (staleness report)
│   ├── cmd_db_setup.py  # `gapit db setup` CLI (bootstrap + materialize + index)
│   ├── mcp.py           # MCP stdio server (hand-rolled JSON-RPC 2.0); backs gapit-mcp
│   ├── mcp_tools.py     # MCP tool implementations (call the shared use-cases)
│   ├── mcp_schemas.py   # MCP tools/list declarations (names, descriptions, inputSchemas)
│   ├── errors.py        # typed errors + JSON error envelope
│   ├── formats/
│   │   ├── tsv.py       # abricate-compatible TSV/CSV (phenotype-free in every typing version)
│   │   ├── csvutil.py   # CSV RFC 4180 + TSV cell sanitization (streaming renderers)
│   │   ├── json.py      # versioned JSON (gapit.report/1 et al.)
│   │   ├── reads_json.py # versioned JSON (gapit.reads/1, gapit.reads/2)
│   │   ├── reads_tsv.py # streaming reads TSV/CSV table
│   │   ├── md.py        # Markdown (human + agent readable, YAML frontmatter)
│   │   ├── cluster.py   # gapit.cluster/1 models + JSON/TSV renderers (typed variant)
│   │   ├── cluster2.py  # gapit.cluster/2 models + JSON/TSV renderers (multi-scheme phenotypes map)
│   │   ├── cluster_md.py # gapit.cluster/1 + /2 Markdown renderers (Phenotype(s) column when typed)
│   │   ├── typing_result.py # gapit.typing_result/1 models + TSV/JSON/MD renderers
│   │   ├── schemas.py   # registered output models behind `gapit schema`
│   │   └── summary.py   # summary matrix renderers (TSV/CSV/JSON/MD)
│   ├── providers/       # 20 database provider modules (13 gene + 7 kaptive cluster) + cluster_common.py + common.py
│   └── py.typed
├── scripts/
│   ├── cluster_calibration.py # developer calibration harness for typed cluster dbs (single-scheme split dbs)
│   └── generate_vp_typing.py # regenerates the committed vpa_ok typing asset from the vpautils OAgc/CPSgc GBKs
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
  `gapit schema report | reads | reads2 | summary | cluster | cluster2 | dblist | dbsetup | dbinfo | dboutdated | datadir | error | version` (plus
  `cluster | cluster2 | features | typing | typing_result | floors`). Keys are snake_case,
  units explicit (`identity_pct`, `coverage_pct`). Semver the schema; never rename or
  retype a field in a minor bump.
- **Markdown**: YAML frontmatter (tool version, db, params, ISO-8601 UTC timestamp) + tables a
  human can read and an agent can regex reliably.
- **Errors**: failures print a JSON envelope to stderr
  `{"schema": "gapit.error/1", "code": "...", "message": "...", "context": {...}}` and exit with a
  documented non-zero code (2 = usage, 3 = missing dependency, 4 = db error, 5 = input error).
- **Typing** `[gapit-extension]`: two-stage designation — `gapit screen -o result.tsv --db NAME`
  detects genes (screen output is phenotype-free; typed and untyped gene dbs screen
  byte-identically), then `gapit typing result.tsv [-D|-f|-o|-q]` designates from the screen
  table (`gapit.typing_result/1`; cluster dbs refuse — their typing rides the screen).
  The MCP `typing` tool accepts the screen table text and returns the same document.
- **DB acquisition** `[gapit-extension]`: `gapit db fetch|path|list|info|search|outdated|build|install|rename|setup` —
  database fetch from the provider catalog (every provider downloads from upstream at fetch
  time — nothing catalog-side is bundled inside the wheel beyond the sanctioned exceptions.
  No provider snapshot rides the wheel, even when its upstream license permits it: CARD/VFDB
  terms are non-commercial, Kaptive is GPL, and provider content carries provenance risk
  beyond the license text. The sanctioned exceptions are the four gapit-curated bundles
  under `src/gapit/data/dbs/` — the public-domain `ecoli_dec`, `lmo_sero`, and `cje_cps`
  panels plus the MIT `vpa_ok` typed cluster snapshot — a pure built-in in the
  ecoli_dec/lmo_sero pattern: no provider module, no `db fetch`; vpautils OAgc+CPSgc merged
  content, MIT (GPL/NC content is never bundled) — materialized into the datadir on first
  use via `bundled.py` (gene bundles replay the gene pipeline, the cluster bundle the
  cluster pipeline); a bundled refresh means a new snapshot riding a wheel release. What
  may ride the wheel is governed by the §6 policy),
  database listing (`gapit db list`, whose TYPING column reads `command` for gene-kind
  databases designated by `gapit typing`, `integrated` for cluster-kind databases whose
  phenotype rides `gapit screen`, and `-` for screen-only databases), per-database `db info -d NAME` (provenance + literature
  references + usage + caveats, merging the provider catalog, `bundled.json`,
  the manifest, provider modules, and an optional `info.json` sidecar;
  `gapit.dbinfo/1`),
  records.jsonl gene search, staleness
  report, custom FASTA→native-db build, GBK/GFF→cluster-db build (`--typing` installs a
  validated `gapit.typing/1` or `/2` spec on either kind), the eight download-on-fetch
  cluster providers (the seven Kaptive panels plus `vch_oagc`, the V. cholerae
  O-antigen/ompW provider; never bundled), the MIT-licensed `vpa_ok` typed cluster bundle (OAgc+CPSgc merged: the committed typing asset
  installs a 3-scheme `gapit.typing/2` doc so one screen renders `gapit.cluster/2` with
  the composed `O:K` serotype, zero network from first install), SHA256-verified
  local-file install, `db build NAME` with **no input** (re-projects an installed gene
  database from its edited `records.jsonl` — the curation path; info/typing/floors
  sidecars and manifest provenance survive, input-shaping flags are rejected), and
  `db rename OLD NEW` (rebuilds an installed database under a new
  name — gene dbs re-header every record, cluster dbs re-stamp the manifest — then removes
  the old directory; the migration path for names that predate the convention).
  **Database naming**: new species databases use the genus initial
  plus the first two letters of the species (`cje` for *C. jejuni*, `eco` for
  *E. coli*, `vch` for *V. cholerae*); genus-level databases use the first three
  letters of the genus (`shi` for *Shigella*, `cam` for *Campylobacter*), followed by a
  content suffix (`cje_cps`, `shi_serotype`). Existing published names are legacy and
  stay stable; `db rename` is the sanctioned migration path when a name must change.
  Datadir precedence is `--datadir`, `$GAPIT_DATADIR`,
  `$XDG_DATA_HOME/gapit/db` (absolute XDG only), then `~/.local/share/gapit/db`;
  `gapit db path [--json]` introspects it without creating, while `db list`/`db setup`
  bootstrap a missing root. See `docs/datadir.md` for the database anatomy, sharing and
  Docker policy.
- **MCP** `[gapit-extension]`: `gapit mcp` / `gapit-mcp` stdio server exposing twelve tools:
  read-only `screen` (incl. `aligner minimap2` assembly survey), `screen_reads` (FASTQ via
  minimap2), `summary`, `typing` (screen-table text → typing result + follow-up notes),
  `schema`, `db_search`, `db_outdated`, `db_path` (resolve the effective datadir without
  creating it), plus `db_list` (may create the empty datadir on first use), and the
  datadir-mutating `db_setup` (bootstrap/index), `db_fetch` (installs databases from the
  provider catalog; may download) and `db_build` (writes a custom db); tool failures carry
  the `gapit.error/1` envelope.
- **stdout purity**: data on stdout, diagnostics on stderr, always. `--quiet` only affects stderr.
- **Self-description**: `gapit --version --json`, `gapit db list --json`, `gapit schema` — an agent
  must be able to discover everything without reading docs.

## 6. Bundled-database policy (what may ride the wheel)

`src/gapit/data/dbs/` is the sanctioned home for license-clean content shipped inside the
wheel: a directory qualifies by carrying `sequences` + `bundled.json` (auto-discovered by
`gapit.bundled`; `kind: cluster` in the metadata additionally requires `features.json`).
Materialization is kind-aware: gene bundles replay the gene-build pipeline, cluster
bundles the cluster-build pipeline (`bundled_cluster.py` rebuilds from the committed
locus FASTA + resolved feature table — gene ids stay snapshot-resolved, the typing spec
re-validates against them).
Bundling stays the exception, not the pattern: every provider fetches from upstream at
`db fetch` time. Third-party provider snapshots are never bundled, even under a
permissive license — the wheel carries gapit-curated content only. Every new bundle
passes three gates, in order:

1. **Curation gate** (hard): the content must be gapit-owned or gapit-curated from
   public-domain sources. A snapshot of a third-party provider database (ncbi, card,
   ecoh, resfinder, ...) does not qualify, regardless of its license.
2. **License gate** (hard, independent of location): only public-domain or permissive
   content (PD, Apache-2.0, BSD, MIT), provenance-audited. GPL (Kaptive) and
   non-commercial terms (CARD, VFDB) never ride the MIT distribution — moving files
   cannot fix this.
3. **Size gate** (budget): the bundled-root total stays ≤ 30 MiB, mechanically enforced by
   `tests/test_bundled.py` (every `pip install gapit` re-downloads all bundles; PyPI caps
   one file at 100 MiB).

Growth tiers past the budget — split rather than inflate the ceiling:

- **≤ 30 MiB**: keep adding to `src/gapit/data/dbs/`. Atomic data-with-code versioning,
  zero packaging machinery; the abricate precedent ships 16 MiB in-package.
- **30–100 MiB**: split a separate `gapit-data` wheel with its own version (the
  `spacy[lookups]` → `spacy-lookups-data` pattern), pulled via a `gapit[datasets]` extra;
  `bundled_databases()` grows a second discovery root.
- **> 100 MiB, or fast-moving upstreams**: fetch-on-demand only. `db fetch` is the
  freshness path for every bundled name anyway — snapshots are point-in-time, and each
  bundled upstream refresh would otherwise force a wheel release.

## 7. Testing strategy

- **Unit**: pure functions (coverage %, merge rules, header parsing) — no I/O beyond `tests/data`.
- **Golden files**: fixed tiny db + fixed contigs → expected TSV/JSON/MD committed; update the
  committed files deliberately and review diffs like code.
- **Parity harness**: `pixi run -e parity parity` (and `summary-parity`) runs real abricate
  (conda) and gapit over a small genome corpus and diffs the gene calls byte-for-byte (file,
  gene, %identity, %coverage). **From v0.5.0 this is a regression reference for the frozen
  baseline surface, not a release gate for new work** — new features need not check abricate
  and may diverge deliberately.
- Full suite + CI matrix Python 3.12; parity byte-diff is opt-in via the `parity` pixi env.
  Offline tests stay fast.

## 8. Domain knowledge (abridged; full detail in SPEC.md)

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
  The default path does **not** merge overlapping intervals (frozen baseline behavior;
  the opt-in `--merge-fragments` cross-contic extension in `fragments.py` is the sanctioned
  exception).
- **Cluster databases** are GBK/GFF-built locus panels screened by minimap2 `asm20` on
  assembly FASTA only. `--r1/--r2` and `--aligner minimap2` are usage errors against a
  cluster database (assemble FASTQ first): locus coverage/identity thresholds are
  assembly-calibrated, and read mapping to repeat-rich multi-kilobase loci needs a
  read-specific preset, depth/MAPQ filtering, and truth-set calibration. Cluster builds
  (and `db setup`) write the `sequences.asm20.mmi` target index; legacy uses it when
  fresh, gseek keeps its own content-addressed cache.
- Key computed fields: `%COVERAGE = 100*(length-gaps)/slen` (filtered unrounded, displayed
  `%.2f`), `%IDENTITY` (BLAST pident, never post-filtered), `COVERAGE_MAP` (15-char minimap),
  `GAPS`. Exact formulas, the dedup rule, and the minimap arithmetic live in `SPEC.md` — they
  are **gapit's frozen contract** for the baseline surface: do not change them silently. New
  surfaces (reads, cluster, typing) follow their own documented contracts and may diverge from
  abricate by design.

## 9. Git & workflow

- Conventional Commits (`feat:`, `fix:`, `test:`, `docs:`, `refactor:`).
- Never commit unless the user explicitly asks. Never force-push.
- `PLAN.md` tracks phases; update it when scope changes, not retroactively.

## 10. For agents working in this repo

1. Read `SPEC.md` before touching `blast.py`, `hits.py`, or `minimap.py`.
2. Public contract changes (JSON schema, exit codes, CLI flags) require a schema version bump and
   a note in `PLAN.md`.
3. Do not add dependencies without checking this file's §2 first — the list is intentionally short.
4. Verification is part of the task: `pixi run lint && pixi run typecheck && pixi run test` before
   declaring done.
