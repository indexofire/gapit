# Changelog

All notable changes to gapit are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `gapit db list` renders a styled rich table (title "Databases"; cyan Provider, status-tinted
  Status, yellow DBTYPE columns) when stdout is an interactive terminal. Piped or redirected
  output keeps the byte-identical TSV, and `--json` is unchanged.
- Provider license metadata surfaced in `gapit db list --json` (`license` per provider when
  pinned) and stamped into the build manifest (card, vfdb, ecoli_vf, kaptive).
- `vendor` field in `gapit.dblist/1` entries: the upstream maintainer organisation behind
  every provider (NCBI, DTU CGE, Kaptive (klebgenomics), ...). Additive; `name` still
  carries the database name passed to `--db`.

### Changed

- Uniform help-on-bare: commands that require input now print their full help when
  invoked with no arguments at all — `gapit screen`, `summary`, `schema`, `db`,
  `db fetch`, `db install`, `db build`, `db search` (typer `no_args_is_help`). The
  help exits 2, the same code the bare root app has always used under click's
  no-args-is-help semantics (recorded reality; the pre-change bare `screen`/`summary`
  also exited 2, but with a `gapit.error/1` envelope instead of help, and bare
  `schema`/`db install`/`db build`/`db search` raised click missing-argument errors).
  Invocations with flags present but input missing keep their typed usage-error
  envelopes (exit 2). Commands valid with no arguments (`mcp`, `setupdb`, `db list`,
  `db outdated`) are unchanged. Breaking for `db fetch`: an omitted NAME no longer
  installs the default set — the default-set download is now explicit
  (`gapit db fetch all`, same DEFAULT_DBS order and receipts; bare fetch prints help,
  and `db fetch --flags` without a NAME is a typed usage error), because a bare
  invocation silently starting a multi-database download was a footgun. The MCP
  `db_fetch` tool aligns: `name` is required and accepts the literal `all`. Docs
  updated in step (EN + zh); no output schemas or goldens changed.

- Terminology sweep after the NAME/PROVIDER column split: user-facing text now reserves
  PROVIDER for the upstream maintainer organisation and consistently calls the `--db`
  value the database NAME. `gapit db fetch` help reads "Fetch and build database(s)"
  with a "Database name" `[NAME]` argument help, `gapit db list` help reads "List known
  databases (NAME, upstream PROVIDER) and their installed state", the `db` group help
  says "database fetch" (was "provider fetch"), `db install` help says "no database
  names" (was "no provider IDs"), and the MCP `db_fetch`/`db_list` tool descriptions use
  database-name phrasing. The unknown-name error from `gapit db fetch` / the `db_fetch`
  tool now reads `unknown database: NAME (available: ...)` with context key `db` (was
  `unknown provider: ...` with context key `provider`); the error CODE stays `USAGE_ERROR`.
  Docs updated in step (EN + zh). No flags, output schemas, exit codes, or internal
  identifiers changed; the `gapit.dblist/1` `providers` array key is untouched.

- `gapit db list` TSV header is now `NAME PROVIDER STATUS DBTYPE DESCRIPTION`: the old
  PROVIDER column actually held the database name users pass to `--db`, so it is renamed
  NAME, and a new PROVIDER column names the upstream maintainer organisation (the `vendor`
  above). Breaking for TSV consumers (header change, 4 → 5 columns); the rich table gained
  the same NAME/PROVIDER split and `gapit.dblist/1` `name` is unchanged.

- card and vfdb unbundled to download-on-fetch providers (license compliance: their
  non-commercial upstream terms — McMaster's for CARD, CC BY-NC for VFDB — are incompatible
  with redistribution inside the MIT-licensed wheel). `ecoli_vf` pins its true license
  (`CC BY-NC 4.0 (VFDB-derived content)` — the repo is labeled Apache but its README states
  the content is taken from the VFDB). The bundled-snapshot archives, the snapshot-first
  install path, the monthly snapshot-refresh CI workflow, and the `snapshot-update` /
  `stale+snapshot-update` `db outdated` statuses are gone; `db outdated` now reports `ok`
  or `stale` only.
- Kaptive cluster providers renamed/expanded to the seven official Kaptive v3 install
  keywords — `kpsc_k`, `kpsc_o`, `kosc_k`, `kosc_o`, `ab_k`, `ab_o`, `ecoli_kps`
  (K. pneumoniae species complex K/O, K. oxytoca species complex K/O, A. baumannii K/OC,
  E. coli group 2+3 capsular polysaccharide loci) — each sourcing its raw GenBank file
  from `main` of the actively curated per-species upstream repos
  (klebgenomics.github.io/Kaptive/db/overview.html#available-databases), so fetches track
  current curation instead of the frozen v2.0.9 archives the four old providers pinned.
  Still GPL-3.0 download-on-fetch, never bundled; manifests record the new source URL,
  license, and citation note (Kaptive/Wyres et al. 2020). The upstream `.toml`
  identity-threshold metadata is not fetched in v1 (a future `typing.json` source).

### Removed

- `gapit db fetch --from-source` flag (meaningless now that no bundled snapshots exist:
  every fetch downloads from upstream).

## [0.4.0] - 2026-09-30

### Added

- **Gene-cluster databases and cluster screening** (stages 1–3; the blastn gene path and
  all existing outputs are untouched — parity stays byte-identical):
  - `gapit db build NAME input.gbk|input.gff3` detects the input kind by suffix and builds
    a `kind: cluster` database — locus FASTA `sequences`, `gapit.features/1` feature table
    (`features.json`), makeblastdb index, manifest `kind: cluster` — with `--kind
    gene|cluster` to override (usage error on contradiction). New parsers (no Biopython):
    GenBank FEATURES/ORIGIN (Bakta/modern and kaptive-style source-note labels) and GFF3
    (embedded `##FASTA` or sidecar); compound `join()` CDS locations are rejected with a
    typed error.
  - `gapit screen --db <cluster db>` dispatches to the minimap2 `asm20` cluster engine:
    per-gene present/partial/absent verdicts (`--min-gene-cov`/`--min-gene-id`, default
    90/90), cs-based union coverage/identity, cross-contig fragmentation resilience, and a
    best-locus call above `--min-cluster-cov` (default 96). Output is the new
    `gapit.cluster/1` document (json), the cluster TSV/CSV (one row per file), or
    Markdown; kind-mismatched flags are usage errors; reads mode rejects cluster
    databases. The MCP `screen` tool gains `minGeneCov`/`minGeneId`/`minClusterCov`.
  - `--typing FILE` installs a validated `gapit.typing/1` declarative scoring spec into a
    cluster db (`typing.json`). Screening a typed database annotates every best call with
    a `phenotype` and an additive `phenotype_detail` breakdown (score, confidence,
    per-rule components, runner-up, ambiguity pair); the TSV gains a PHENOTYPE column
    (typed runs only — untyped output stays byte-identical). Three rule kinds —
    `weighted_genes`, `cluster_match`, `learned_linear` (schema only; no training code) —
    behind a `cutoff` + `ambiguity_margin` + `fallback` decision layer. Unknown gene/locus
    references fail with `TYPING_UNKNOWN_GENE` at build AND screen time.
  - Four kaptive cluster providers, `kaptive_k`, `kaptive_o`, `kaptive_ak`, `kaptive_oc`
    (Klebsiella K/O and A. baumannii K/OC antigen loci; Wyres et al. 2020 — cite Kaptive).
    The database content is GPL-3.0, so nothing is bundled: fetch downloads the pinned
    Kaptive v2.0.9 GenBank files and builds cluster databases; manifests record the
    source URL, `license`, and citation `note`. No typing model ships with them yet
    (locus calls only, phenotype null).
  - `scripts/cluster_calibration.py`: developer harness that screens a typed cluster db
    against a labels TSV and prints the per-expected-phenotype score distribution, the
    called×expected agreement matrix, and the divergence list (the v1 loop for future
    learned_linear training).
- `gapit screen --merge-fragments` (blastn contig mode only; off by default): merge gene
  fragments split across contig boundaries into one reported hit when their union subject
  coverage reaches `--mincov`. Merged rows carry additive optional `merged` and `fragments`
  fields in `gapit.report/1` JSON plus a Markdown detail line; the default path (and all
  existing output) is untouched. The MCP `screen` tool gains the matching `mergeFragments`
  boolean parameter.

## [0.3.1] - 2026-09-28

### Added

- `-h` now works as an alias for `--help` on every command (root, subcommands,
  and nested `db` sub-app).

## [0.3.0] - 2026-09-28

### Removed

- `gapit list` command and the `gapit.list/1` schema (breaking; rightsholder decision to
  drop abricate `--list` strict parity for this surface). `gapit db list` / `gapit.dblist/1`
  is the single listing surface; `gapit schema` now introspects six documents.

## [0.2.2] - 2026-09-24

### Security

- minimap2 input paths are passed as absolute paths, so a file named like an option
  (e.g. `-d`) can no longer inject minimap2 flags; the child process no longer inherits
  stdin (the MCP protocol stream).

### Changed

- `--minid`/`--mincov` with `--aligner minimap2` or `--r1/--r2` are now a usage error
  (exit 2) instead of being silently ignored.
- MCP `screen` (blastn) reuses the CLI use-case; its validation messages now match the CLI.
- MCP `screen_reads` passes read arrays natively, so filenames containing commas work.

### Fixed

- MCP server replies `-32600`/`-32602` to malformed requests that carry an `id` instead of
  silently dropping them.
- A truncated gzip reads file raises a typed `INVALID_READS_FORMAT` error (exit 5) instead
  of `UNEXPECTED` (exit 1).

## [0.2.1] - 2026-09-22

### Changed

- License changed from GPL-2.0-only to MIT (rightsholder decision; behavioral
  reimplementation status unchanged).

## [0.2.0] - 2026-09-22

Development window: 2026-09-15 to 2026-09-22 (PLAN.md Phases 1 through 8 and post-phase
hardening).

### Added

- Screening core: `any2fasta` to `blastn` pipeline with SPEC-exact hit processing, and
  byte-compatible TSV against abricate 1.4.0 (PLAN Phases 1-3, 2026-09-15). Parity is
  checked on a genome corpus with `pixi run -e parity parity` (6/6 byte-identical).
- Agent-facing outputs: `--format json|md` with versioned schemas (`gapit.report/1`,
  `gapit.reads/1`, `gapit.summary/1`, `gapit.list/1`, `gapit.version/1`), the
  `gapit.error/1` stderr envelope with documented exit codes, and introspection via
  `gapit schema report|reads|summary|list|error|version` (2026-09-16).
- FASTQ reads mode via minimap2: `gapit screen --r1/--r2 --read-type sr|map-ont|map-hifi`;
  presence is called on alignment breadth (2026-09-16).
- Summary matrix mode: `gapit summary` emits TSV/CSV/JSON/MD; byte parity with
  `abricate --summary` via `pixi run -e parity summary-parity` (2026-09-17).
- `gapit db install`: checksum-verified (SHA256), atomic local-file installation
  (2026-09-17).
- Native `gapit/v1` database format: tagged-header codec, `records.jsonl` truth store
  with `gapit.manifest/1` provenance, deterministic build pipeline with self-check, and
  `.mmi` BLAST index build with version-gated reuse (2026-09-18).
- 12 database providers (ncbi, card, resfinder, argannot, plasmidfinder, megares, ecoh,
  vfdb, ecoli_vf, bacmet2, victors, upec_expec_vf) with offline-tested transforms
  (2026-09-18).
- Bundled card and vfdb snapshots: a bare `gapit db fetch` installs the default set with
  zero network; `--from-source` forces upstream download (2026-09-18).
- CLI hardening: shell completions, `--debug` argv echo on stderr, and `--jobs` parallel
  screening across inputs with input-order stdout (2026-09-18).
- CI: gates matrix (lint, fmt, typecheck, test on Python 3.11 / 3.13 / 3.14) plus a
  parity job that byte-diffs against real abricate (2026-09-18).
- MCP stdio server: `gapit mcp` subcommand and `gapit-mcp` console script, hand-rolled
  JSON-RPC 2.0 with zero new dependencies; read-only tools `screen`, `summary`,
  `schema`, `db_list` (2026-09-18).
- Conda recipe (`recipe/meta.yaml`) for PyPI/bioconda packaging (2026-09-18).
- Custom database construction: `gapit db build NAME FASTA` builds a screening-ready native
  database from plain, abricate `~~~`, or `gapit|` FASTA (auto-detected), with optional
  `--tsv` metadata (accession, function classes) and `--dbtype` override; guide with worked
  examples in `docs/custom-db.md` (2026-09-19).
- `--aligner blastn|minimap2` engine selector on `gapit screen`: defaults follow the input
  (blastn for contig files, minimap2 for `--r1`/`--r2` reads); `--aligner minimap2` screens
  positional assembly FASTA through the minimap2 engine, `--aligner blastn` with `--r1`/`--r2`
  is a usage error (2026-09-19).
- `gapit db outdated`: staleness report over installed databases (age vs `--days`, newer
  bundled snapshot) as a TSV table or `gapit.dboutdated/1` JSON document (2026-09-20).
- `gapit db search TERM`: case-insensitive gene/accession/function/product lookup across
  installed databases' `records.jsonl`, with `--field`, `--exact`, `--db`, `--limit` and
  JSONL output (2026-09-20).
- MCP database tools `db_fetch`, `db_build`, `db_search`, `db_outdated`: the MCP server now
  exposes the `db` commands alongside the analysis tools, so agents can self-provision and
  inspect databases mid-session; `mcp.py` split into protocol (`mcp.py`) and tool
  implementations (`mcp_tools.py`) with zero wire change for the existing tools (2026-09-20).
- CI snapshot refresh: monthly scheduled workflow (`.github/workflows/snapshot-refresh.yml`)
  re-fetching card and vfdb from upstream, regenerating the bundled tars only when the
  records changed, and opening a review PR (2026-09-20).
- Reads-mode opt-in alignment filtering: `gapit screen --min-identity/--min-mapq` (reads
  engine only) drops PAF alignments below the thresholds before aggregation and emits the new
  `gapit.reads/2` document (params gain both thresholds; gene entries gain
  `mean_identity_pct`, the alignment-length-weighted mean identity; introspectable via
  `gapit schema reads2`). Both flags off keeps `gapit.reads/1` byte-identical, and the
  identity floor fixes the documented family-splitting over-calls on homologous genes
  (2026-09-20).
- MCP reads tools: new `screen_reads` tool (FASTQ lanes via minimap2 → `gapit.reads/1`,
  `min_identity`/`min_mapq` > 0 → `gapit.reads/2`) and new `screen` arguments
  `aligner`/`min_breadth`/`min_identity`/`min_mapq` (`aligner minimap2` = fast assembly
  survey). Both delegate to the shared reads use-cases with `quiet` stderr; additive
  `gapit.mcp` contract — nine tools, the other eight wire-identical, no output-schema
  version bump (2026-09-22).

### Changed

- `any2fasta` is no longer a runtime dependency: input normalization (FASTA/FASTQ/GenBank/EMBL,
  plain/gz/bz2 → FASTA) is native (`seqconvert.py`, perl-extracted semantics) and blastn reads the
  converted FASTA on stdin. Parsed-record equivalence with the binary is differentially tested
  with the parity env on PATH (abricate provides any2fasta transitively), and `.fa` parity
  remains byte-identical (2026-09-21).
- Screening and reads paths route `gapit/v1` tagged headers to the native codec; legacy
  abricate `~~~` headers keep the frozen parser, so abricate-built datadirs screen
  identically. The reverse does not hold: gapit-built databases are unreadable by
  abricate (accepted trade-off).
- A bare `gapit db fetch` installs the bundled default set (card, vfdb) instead of
  downloading; `--from-source` restores the upstream fetch path.
- `gapit setupdb` honors the `gapit.manifest/1` dbtype on reindex (manifest-less abricate
  dirs keep the mol_type heuristic), and reads/minimap2-assembly screening now rejects the
  blastn-only flags `--fofn`, `--noheader`, `--nopath`, `--jobs N>1` with usage errors
  (exit 2) instead of silently ignoring them (2026-09-22).

### Removed

- `--csv` flag on `gapit screen`; use `--format csv` instead (breaking: the flag is now
  rejected as an unknown option, exit 2) (2026-09-19).

### Fixed

- Summary mode auto-detects TSV vs CSV per input file, avoiding the upstream quirk where
  a global `--csv` flag mangles mixed-format input sets.
