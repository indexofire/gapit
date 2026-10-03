# Changelog

All notable changes to gapit are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.1] - 2026-10-03

### Added

- **Illumina FASTQ wildcard support: `gapit screen -d ecoli_dec *.gz` screens a shell glob
  of reads without `--r1`/`--r2`** (rightsholder batch workflow). When every positional file
  is FASTQ — a `.fastq`/`.fq` (± `.gz`) extension, or content sniffed as FASTQ when the
  extension is ambiguous (via the existing `detect_read_kind`; known contig extensions and
  garbage content keep today's contig pipeline untouched) — the command enters reads mode
  and auto-pairs samples from filenames (`gapit.readpairs`, new): mate markers `_R1_001`/
  `_R2_001` (bcl2fastq; an `_L00x` lane tag stays in the sample key), `_R1`/`_R2`, `_1`/`_2`,
  `.1`/`.2`, all case-insensitive, longest suffix first, extensions `.fastq`/`.fq` ×
  `.gz`/`.bz2` stripped before matching. The document carries one `files[]` entry per sample
  (`reads[]` holds the sample key; per-gene metrics union over the sample's lanes; samples
  sorted lexicographically). A file whose mate is missing, or with no detectable marker,
  screens single-end with a quiet-respecting stderr `WARNING: no mate found for X —
  screening single-end` — a warning, never an error. Duplicate identical basenames merge as
  extra lanes of one sample. Mixing FASTA and FASTQ positionals is a usage error naming the
  reads files; the reads-mode guards fire as on `--r1`/`--r2` (`--minid`/`--mincov`,
  `--merge-fragments`, ... rejected), while `--min-identity`/`--min-mapq` now work on this
  path (gapit.reads/2). `--read-type` keeps
  the content-based preset resolution (`sr` for FASTQ). **`--jobs` is legal on this path
  only** (the `--r1`/`--r2` single-sample invocation keeps the frozen rejection): samples
  screen concurrently through a ThreadPoolExecutor whose positional yields keep output in
  sample order byte-identical to `-j 1` (a failing sample raises at its position; the
  `--jobs × --threads` oversubscription stderr note and the `--jobs >= 1` validation mirror
   the contig path). **Output defaults to streaming tsv on this path**: the `#SAMPLE` header
   leads, then each sample's gene rows land the moment its screening completes (`--format
   csv` is the comma spelling), so a terminal shows per-sample results as they arrive; json
   stays the opt-in single document written at the end (`--format json`; since the same
   cycle's format-unification it is the opt-in on the `--r1`/`--r2` path too — see Changed).
   md emits
  a static frontmatter first (no run totals on this path — totals live in the JSON
  document; the `--r1`/`--r2` render keeps its `files:`/`genes_found:` frontmatter
  unchanged) then one `## <sample>` section per completed sample, head-of-line in sample
  order under `--jobs`; `--output` writes the streamed chunks incrementally, so a
  mid-batch failure persists samples 1..k-1 (tsv/md) and nothing at all (json). The
  single-sample `--r1`/`--r2` path, the blastn contig path, and
  the explicit `--aligner` routings are unchanged.
- **`gapit typing` closes every TSV/MD row with a `GENES` column and `gapit.typing_result/1`
  gains the additive per-file `genes` list** (schema stays `/1`; additive field, the
  0.5.0 shape unchanged otherwise): the FILE's present gene names — the folded best hits the
  designation ran on — sorted alphabetically and `;`-joined, repeated on each of the FILE's
  scheme rows so each call is auditable against the hits that drove it. JSON carries the same
  list as `files[].genes`.
- **`gapit summary` reads its report table from stdin** — the pipe twin of `gapit typing`:
  `gapit screen -d ecoli_dec *.fna --quiet | gapit summary` summarizes the batched screen
  output straight off the pipe (a piped table is one input, so dutch mode applies — one row
  per FILE value, byte-identical to summarizing the `-o`-written table). Stdin is read
  whenever stdin is not a terminal; a bare invocation at a terminal keeps printing help;
  `gapit summary -` is the explicit stdin marker (reads until EOF even under a terminal);
  mixing `-` with file arguments is a usage error; non-UTF-8 stdin is the `SUMMARY_MALFORMED`
  envelope naming `-`. The flag-but-no-file invocation (`summary --quiet` with a pipe) now
  reads stdin instead of the former USAGE_ERROR.

- **`--db` now also answers to the word-style alias `-db`** (on `gapit screen` and
  `gapit db search`): `-db ecoli_dec` parses exactly like `--db ecoli_dec`, matching the
  tar-style single-dash habit that previously misparsed as `-d b` + a positional. The
  space-separated form is the supported one — click binds an attached `-dbecoli_dec` to
  `-d` with value `becoli_dec` (the original failure), which the new guard below catches.

### Fixed

- **A positional argument that names a known database but is missing on disk now fails with
  a targeted usage error (exit 2)** instead of the generic INPUT_NOT_FOUND (exit 5) plus a
  starved downstream pipe: `gapit screen -d b ecoli_dec X.fna` (the stripped `-db` misparse
  shape) reports `input file not found: 'ecoli_dec' — it is a database NAME; screen takes
  databases via --db/-d/-db (e.g. --db ecoli_dec), positional arguments are genome files`.
  Database names are matched against the provider catalog, bundled wheel content, and
  installed datadir databases; ordinary missing files keep today's INPUT_NOT_FOUND.

### Changed

- **Breaking: `gapit screen` reads mode now defaults to streaming tsv on EVERY path** —
  `--r1`/`--r2` and `--aligner minimap2` join the positional wildcard (and the blastn contig
  path) with tsv as the human default; json/md are the explicit agent opt-ins. tsv is the
  human default on every surface; json (`gapit.reads/1`/`2`) and md are the opt-ins for
  agents. The default output is the same streaming `#SAMPLE` table the wildcard emits (one
  chunk for the single sample; preamble + chunk compose exactly the buffered render);
  `--format csv` is the comma spelling, `--format json`/`--format md` unchanged. The former
  reads json default shipped ≤0.5.0 — invocations that relied on it must now pass
  `--format json` (the tsv/csv reads-mode usage rejection is gone with it). The MCP
  `screen`/`screen_reads` tools are unchanged: they pass `format` explicitly and still
  return json by default (the agents' surface).

### Added

- `gapit summary --coverage` / `-c`: cells show %COVERAGE values (`;`-joined), reproducing the
  classic abricate `--summary` cell shape byte-for-byte (the parity harness now runs in this
  mode). `-ic` combines both metrics: each hit renders `identity/coverage`. `gapit typing` TSV/MD: the `GENES` column moved to follow `PHENOTYPE`
  (was last; additive JSON `files[].genes` unchanged).

### Changed

- **Breaking**: `gapit summary` default cells are now the presence call (`+` / `-`) instead of
  %COVERAGE values — pass `--coverage` to restore the previous cell shape.

### Changed

- **Project stance**: from v0.5.0 onward gapit pursues its own contract — abricate parity is no
  longer a design goal. The contig TSV stays a frozen historical baseline (kept for existing
  pipelines); the parity harness is a regression reference, not a release gate. New surfaces may
  diverge from abricate by design. (AGENTS.md updated accordingly.)
- Reads tables: the default lists present genes only (`--all-genes` / `-A` adds absent
  calls) and the last column is `PRODUCT` (the `RESISTANCE` column is dropped from reads
  tables; the contig table keeps it).

### Fixed

- **Reads-mode breadth accuracy (26ECO0071 astA)**: minimap2 now always runs with `--cs`
  and the `sr` preset carries `-k11 -w6 -B2`. Two stacked defects silently dropped
  diverged-allele genes from FASTQ screening (the assembly/BLAST path found astA, reads did
  not): (1) the stock `sr` seed (k=21) finds no exact run on a ~90%-identity 117 bp gene —
  zero alignments regardless of scoring; (2) minimap2's no-CIGAR mode emits clipped
  alignment coordinates (57% vs 100% breadth on the same data). Presence calls and typing
  inputs derived from reads screening change accordingly; the reads goldens were regenerated
  deliberately.

## [0.5.0] - 2026-10-02

### Breaking

- **`--format md` now streams per file, and the Markdown layout lost its run totals (MD goldens
  and docs regenerated).** The Markdown frontmatter was dynamic — `files:`/`hits:` counts (and
  the cluster report's all-files summary table) only exist once every file has screened, which
  forced the renderer to buffer the whole document until the end. The frontmatter is now STATIC
  metadata only (`schema`, `tool`, `created_at`, `db`, thresholds — everything knowable before
  file 1): it prints first and each file's section (heading, table, merge-fragment lines; on
  cluster databases the file's own summary row plus its gene table, replacing the global summary
  table) prints the moment that file completes — head-of-line in input order under `--jobs`,
  exactly like the tsv rows, with `--output` flushing per chunk. Run totals now live only in the
  JSON document (`gapit.report/1` `files[]`/`hits`, `gapit.cluster/1`), which stays a single
  buffered document by design. Update any parser that read `files:`/`hits:` from Markdown
  frontmatter or expected the cluster summary table above the sections; per-file content is
  unchanged byte-for-byte. Reads-mode Markdown (one document per run by construction) and
  `gapit typing`'s Markdown are unchanged.

- **`gapit screen --db` is now required (no default).** The silent `ncbi` default could
  auto-materialize a bundled 8.4k-record database as a hidden side effect of a bare
  `gapit screen contigs.fa`; a screen now refuses to run (typer usage error, exit 2) until
  `--db NAME` is given explicitly. The MCP `screen`/`screen_reads` tools mirror the change:
  `db` moved into each inputSchema's `required` array (default removed), and an absent
  `db` returns the `gapit.error/1` USAGE_ERROR envelope. Update any invocation that relied
  on the default.
- **`gapit screen` no longer embeds phenotypes — designation moved to the new `gapit typing`
  command (one-cycle revert; the embedded surfaces never released).** Screening a typed gene
  database is now pure gene detection, byte-identical to screening an untyped one in every
  format: the earlier-this-cycle `PHENOTYPE` column on the report TSV/CSV, the
  `files[].phenotypes` object of `gapit.report/1`, and the Markdown phenotype lines are all
  gone (the typing engine itself — rules, schemes, gates, compose — is unchanged, and
  `typing.json` still ships with databases). The reads-mode guard that rejected typed gene
  dbs (`gene-path typing is contig-mode only`) is gone with it: `--r1/--r2` and
  `--aligner minimap2` screen typed gene databases like any other. The cluster path is
  untouched — typed cluster databases keep their integrated designation (the
  `gapit.cluster/1` phenotype block and the typed-only TSV `PHENOTYPE` column). The MCP
  `screen` tool returns the pure report on typed gene dbs (a typing MCP tool is future
  work). If you consumed `files[].phenotypes` from pre-release builds, switch to
  `gapit typing` (below).

### Added

- **`gapit typing RESULT.tsv [RESULT2.tsv ...]` — designation from screen results (the
  rightsholder's two-stage CLI design).** Stage 1: `gapit screen -o result.tsv --db NAME`
  detects genes; stage 2: `gapit typing result.tsv` reads the gapit/abricate screen table
  (TSV or auto-detected CSV, several tables merge by their FILE column), resolves NAME from
  the rows' shared `DATABASE` column under `--datadir` (bundled databases still materialize
  on first use), folds every `(FILE, GENE)` to its best row by `(%IDENTITY, %COVERAGE)`
  (first row wins ties — the same per-gene fold the report hits received), and runs the
  existing typing engine over each FILE's genes. Output: default TSV with columns `FILE`,
  `SCHEME`, `PHENOTYPE`, `CONFIDENCE`, `SCORE`, `RUNNER_UP`, `NOTES` (ambiguous calls
  render `-` and carry the candidate pair in NOTES; streamed per input file), `--format
  json` as the new versioned **`gapit.typing_result/1`** document (`{schema, tool,
  created_at, source, db, files[].phenotypes}` — full score breakdowns, introspectable via
  `gapit schema typing_result`, also registered in the MCP `schema` tool), and
  `--format md` (frontmatter + the same seven-column table); `-o FILE`/`-q` behave as
  everywhere. Typed errors: mixed `DATABASE` values `DATABASE_MISMATCH`, a table with no
  data rows `TYPING_NO_DATA` (a screened-no-hits file has no rows — the database is only
  knowable from them), a database without a `typing.json` `TYPING_NO_SCHEME`, and a cluster
  database `TYPING_CLUSTER_DB` (its typing is integrated in `gapit screen`). Equivalence is
  locked end to end: the DEC matrix, serotyping fixtures, and the bundled `ecoli_dec`
  assert the command's calls equal the calls the inline engine produced before the move,
  with new goldens (`typing_markers.*`, `typing_serotyping.tsv`, `ecoli_dec_typing.tsv`).
  The table also arrives on stdin: `gapit screen 1.fna --db NAME | gapit typing` types the
  piped output with no file arguments (stdin must not be a terminal — a bare invocation at
  a terminal still prints help), the explicit `gapit typing -` marker reads stdin even
  under a terminal (until EOF), mixing `-` with file arguments is a usage error (v1:
  either stdin or files), a piped table's JSON `source` renders as `["-"]`, and an empty
  pipe is the usual `TYPING_NO_DATA` naming `-`. The three delivery forms (file, `-`,
  pipe) are locked byte-identical.
- **Short single-dash aliases for every CLI option (rightsholder UX directive).** Each option
  on every command now pairs a short form with its long form — `--help`/`-h` style, e.g.
  `gapit screen contigs.fa -d ncbi -f json` for `--db ncbi --format json`. Letters follow the
  first-letter convention with uppercase on collision (`--db`/`-d` vs `--datadir`/`-D`;
  `--minid`/`-i` vs `--min-identity`/`-I`), mnemonics where the first letter is taken
  (`--min-gene-cov`/`-g`, `--read-type`/`-x`, `--r1`/`-1`), and `-v` for `--debug` (verbose).
  Deliberate exceptions: boolean negative halves stay long-only (`--no-merge-fragments`) and
  the typer-managed completion flags are untouched. Long forms are unchanged and remain the
  documented canonical spelling; a regression test walks the built command tree proving every
  option has a unique per-command short and byte-compares short vs long invocations.
- **`gapit screen --output PATH` writes the report to a file (all engines, including
  reads mode).** The file opens on the first output byte — truncating any existing file
  (v1 overwrite semantics, never append) — every streamed chunk is flushed, and stdout
  then carries no data (stderr diagnostics unchanged). An in-batch failure on file *k*
  leaves files 1..*k-1*'s already-streamed output persisted in the file (or printed, on
  stdout) before the typed error envelope and documented exit code; a run failing before
  any output (usage/dependency/db errors) creates no file.
- **Per-file streaming emission for tsv/csv.** `gapit screen` (blastn gene path and
  cluster path) now emits results in input order as each file completes instead of
  buffering everything until the end: the TSV/CSV header prints once screening starts and
  each file's rows print the moment its file finishes — head-of-line under `--jobs N > 1`
  (file *i* waits for 1..*i*, so output bytes stay identical to sequential runs).
  Markdown sections are computed per file but the document (frontmatter totals) is emitted
  once at the end; JSON stays a single buffered document. Chunk concatenation is
  byte-identical to the former buffered render — goldens and abricate parity are
  untouched. The use-cases gained an optional `emit: Callable[[str], None]` sink (default
  `None` = buffered return exactly as before, so MCP and library callers are unchanged).

- **Bundled databases — four audited provider snapshots join `ecoli_dec` in the wheel
  (five bundles total).** `ncbi` (8373 records, public domain), `resfinder` (3206,
  Apache-2.0), `ecoh` (597, BSD-3-Clause), and `upec_expec_vf` (77, MIT) now ship under
  `src/gapit/data/dbs/<name>/` as point-in-time snapshots produced by the `db fetch`
  provider pipeline itself (gapit/v1 headers; snapshot date recorded in `bundled.json`
  as `snapshotted`, the metadata model's new field). `db fetch <name>` remains the
  fresh-upstream update path for all four: it re-downloads the latest upstream content
  and (with `--force`) overwrites the materialized datadir copy, whose manifest then
  drops the `bundled` source stamp; a regression test pins the snapshot↔fetch
  byte-identity for `resfinder` against a committed copy of the exact upstream archive,
  and the provenance guard now scans ALL five bundles' headers for `VF*`/`VFDB`/`ARO:`
  markers. `db list` merges registry-and-bundled names into ONE bundled-section row
  (STATUS `bundled` → `installed (N)`; JSON `source: "bundled"`); registry rows for
  non-bundled providers stay byte-identical. Wheel size 170 KiB → 1.5 MiB (+1.3 MiB for
  the four snapshots; their `sequences` are ~13 MiB uncompressed). `db outdated` stays
  manifest-based (bundled-vs-upstream staleness comparison is future work).
- **Bundled databases — `ecoli_dec` ships install-time ready (zero-network first use).**
  A new mechanism (`gapit.bundled`) lets license-clean, public-domain databases ride the
  wheel under `src/gapit/data/dbs/<name>/` (`sequences` + optional `typing.json` +
  `bundled.json` metadata: name, description, vendor, dbtype). The first database to use
  it is `ecoli_dec`: the 17-record diarrheagenic E. coli marker panel with the dual-scheme
  `gapit.typing/2` designation (gb4789_6 + risk_monitoring) — after the 2026-10 provenance
  audit re-sourced all records to NCBI/DDBJ primary submissions. When `--db` names a
  bundled database that is absent from the datadir, the screen path materializes it first
  (one quiet-respecting stderr note, then the standard gene-build pipeline: records.jsonl,
  gapit/v1 sequences, BLAST index, typing copy, manifest stamped with the additive
  `source: "bundled"` field); a completely missing datadir is bootstrapped on this path
  only. `gapit setupdb` materializes every bundled database alongside indexing; both paths
  are idempotent (no-op once the manifest exists). `db list` gains a bundled section
  between the registry rows and local extras on all four surfaces (TSV, rich table,
  `--json`, MCP `db_list`): STATUS `bundled` before materialization, `installed (N)`
  after, JSON `source: "bundled"` — registry and local rows stay byte-identical. Wheel
  size 159 KiB → 170 KiB (+11 KiB for the panel).
- **DEC dual-scheme designation + the `requires_any` exact_set primitive.** The
  diarrheagenic E. coli placeholder is replaced by the real designation as TWO schemes over
  one gene database (`tests/data/typing/schemes/dec.json`, the typing/2 multi-scheme
  showcase): `gb4789_6` (GB 4789.6-2016 panel semantics — EAEC is any-of aggR/pic/astA)
  and `risk_monitoring` (最新食品安全风险监测方案 — aggR mandatory), both exact
  semantics (cutoff 1.0, margin 0.0) with the `uidA` control gene, the `non-DEC` fallback,
  and one severity-ordered rule ladder (EHEC > STEC/EPEC > ETEC > EIEC > EAEC; hybrids
  surface as the runner_up; `EPEC_atypical` declared after EHEC/STEC so stx+ isolates never
  land there). To encode it, `exact_set` rules gained an optional `requires_any` any-of
  gate (satisfied iff all `requires` present AND, when the set is non-empty, at least one
  `requires_any` present AND no `excludes` present; `requires` is now optional when
  `requires_any` carries the constraint; the three sets must be pairwise disjoint and not
  all empty — `TYPING_MALFORMED` otherwise; additive within typing/2, no schema bump).
  The synthetic fixture `dec.fa` carries all 14 gene names plus duplicate pic/sth records
  mirroring the rightsholder db (duplicate records collapse by gene name:
  `-culling_limit 1` keeps the best subject, the engine folds one call per name). The
  definitional matrix is locked end to end — including the headline divergence
  (pic+astA+uidA → EAEC under gb4789_6, non-DEC under risk_monitoring, golden
  `gene_dec_pic_astA`) — and hybrid severity ties (stx2a+escV+aggR → EHEC with EAEC at
  1.0 as runner_up). No output-schema, TSV, or parity changes.
-   **typing/2 cycle complete — the scheme cookbook (stage 3)**, consolidating the typing/2
  story: stage 1 promoted the typing document to `gapit.typing/2` named multi-scheme
  documents (a `/1` document degrades to one anonymous `default` scheme; gene-kind builds
  gained `--typing`, one call per scheme — first keyed into gapit.report/1's phenotypes
  object, later moved to the `gapit typing` command), stage 2 added the six rule/scheme
  primitives (entry below), and stage 3 proves the framework encodes real designation
  schemes: six researched schemes as validated example documents with synthetic marker
  fixtures in `tests/data/typing/schemes/`, each self-checked end to end by the new
  `test_typing_schemes*.py` suites (build fixture db → screen → asserted phenotype calls)
  with goldens for one representative call per scheme (doumith 4b, shigella mixed, the
  vp O:K compose, cholerae inaba). The cookbook: the complete Doumith Listeria table
  (Doumith et al. 2004, incl. the 4b*/IVb-v1 HGT caveat falling back to NT), a
  ShigaTyper-semantics Shigella/EIEC skeleton (Wu et al. 2019: ipaH_c control, wzx
  unique-group mixed calls, sonnei form I/II with the exact-tie ordering, the lacY EIEC
  approximation with the S. boydii 9/15 exemptions), the meningotype serogroup panel
  (Mothershed et al. 2004 + the synG EX7E allele-probe trick at identity_floor 99.5),
  the Vibrio parahaemolyticus O/K Kaptive pattern (van der Graaf-van Bloois et al. 2023:
  cluster_match bands, the O3/O13 combined label, OUT/KUT fallbacks, per-scheme databases
  — a cluster db carries exactly one scheme), V. cholerae O1/O139 + Ogawa/Inaba (wbeT
  single-SNP allele floor, the negative-wbeT Inaba rule, Hikojima not determinable, the
  wbfZ junction-gene trap note), and a documented DEC placeholder awaiting the
  rightsholder's curated panel. New bilingual docs page `docs/typing.md` / `typing.zh.md`
  (framework reference, primitive tables, the cookbook with citations, the allele
  probe-trick section and its future allele_match successor, calibration pointer), wired
  into mkdocs nav (EN + zh), the doc indexes, and the databases pages. No schema, output,
  or TSV parity changes; parity 15/15 + summary 6/6 unchanged.
- typing/2 stage-2 rule primitives for declarative phenotype scoring (additive; existing
  rule semantics, outputs, and goldens untouched): the `exact_set` rule (deterministic
  boolean gene-set match with `requires`/`excludes` and 90/90-default identity/coverage
  floors — a 1.0 tie involving a satisfied exact_set resolves by declaration order
  instead of the ambiguity margin), optional `coverage_floor` on `weighted_genes` and
  `exact_set`, scheme-level `control_gene` (absent control gene zeroes the scheme to its
  fallback with a "control gene absent" note), scheme-level `unique_group` +
  `mixed_phenotype` (multiple present group genes call the mixed phenotype with the pair
  in `ambiguous`), `compose` schemes rendering `"{o_group}:{k_group}"` from sibling
  schemes' calls (fallback strings flow through; an ambiguous ingredient yields a null
  composition carrying that ingredient's variants), and rule-level `notes` surfaced
  verbatim on winning phenotype calls. All new fields validate (`TYPING_MALFORMED`,
  unknown gene references `TYPING_UNKNOWN_GENE`); rule models split into
  `gapit.typing_rules` and the decision layer into `gapit.typing_decide` at the 250-LOC
  ceiling.
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

### Fixed

- **DEC panel provenance cleaning.** The original 17-record DEC panel carried two
  VFDB-derived records (spotted past header inspection) plus VF-flavored metadata; all
  records were re-sourced from NCBI/DDBJ primary submissions (100% same-allele, 5
  replacements: bfpB, pic ×2, astA, stx2b, and the sth/stp/lt accessions re-verified),
  and a regression test now mechanically asserts no `VF*` tags ride the bundled headers.

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

- `gapit db list` now also lists databases installed into the datadir outside the provider
  catalog: `db build` products (gene and cluster kinds) and manifest-less directories
  (abricate-style or `db install` bytes) render after the registry entries, sorted by name,
  with PROVIDER `local`, the record count from the manifest (FASTA count when absent), and
  DBTYPE from the manifest / BLAST index suffix / letter heuristic. The gapit.dblist/1 JSON
  gains them inside `providers` behind the additive `source: "local"` field — registry
  entries simply lack the field, so pre-existing output is byte-identical; TSV, rich table,
  and the MCP `db_list` tool share the same rows.
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
