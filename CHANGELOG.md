# Changelog

All notable changes to gapit are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
