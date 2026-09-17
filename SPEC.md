# SPEC.md — abricate behavior specification (parity contract for gapit)

> Distilled from abricate **1.4.0**, master commit
> [`e2064df7d193ad783d4c188d4ce79706faf9eb75`](https://github.com/tseemann/abricate/commit/e2064df7d193ad783d4c188d4ce79706faf9eb75)
> (2026-07-17). abricate is a single 530-line Perl script (`bin/abricate`) plus a DB builder
> (`bin/abricate-get_db`); there are no perl5 libraries.
>
> **Parity principle**: where gapit and this document disagree with intuition, abricate wins.
> Every rule below is a parity requirement unless marked `[gapit-extension]`.

## 1. CLI surface

| Option | Type | Default | Notes |
|---|---|---|---|
| `--db` | str | `ncbi` | subdir of datadir |
| `--datadir` | path | `<script-dir>/../db` | gapit: default to env var `GAPIT_DATADIR`, then platform data dir |
| `--minid` | float | **80** | `0 < minid <= 100`; enforced only via blastn `-perc_identity` |
| `--mincov` | float | **80** | `0 <= mincov <= 100`; post-filter on unrounded float |
| `--threads` | int | 1 | passed to `-num_threads` |
| `--fofn` | path | — | file-of-filenames; **replaces** positional args |
| `--quiet` | flag | off | silences stderr only |
| `--csv` | flag | off | field separator `,` instead of tab |
| `--noheader` | flag | off | suppress `#FILE ...` header row |
| `--nopath` | flag | off | basename the FILE column |
| `--summary` | flag | off | summary-matrix mode (§7), args are report files |
| `--identity` | flag | off | summary cells show %IDENTITY instead of %COVERAGE |
| `--list` / `--setupdb` / `--check` / `--version` / `--help` | modes | — | |
| `--debug` | flag | off | verbose stderr |
| `--format` `[gapit-extension]` | enum | `tsv` | `tsv\|csv\|json\|md`; supersedes upstream's validated-but-unimplemented `--outfmt` (`bed gff json` are accepted upstream but do nothing) |

Mode precedence (upstream): `--summary` → `--check` → dep check → `--list`/`--setupdb` → BLAST
version gate (`blastn -version` must be ≥ 2.2.30; we require modern BLAST+ ≥ 2.7 via conda) → run.

Exit codes (upstream): `0` ok; `1` any runtime error; `5` unknown option; BLAST/any2fasta pipeline
failure propagated verbatim. **gapit mapping** `[gapit-extension]`: `2` usage, `3` missing
dependency, `4` db error, `5` input error, `1` unexpected — stdout TSV stays byte-compatible;
exit-code integers are gapit's own contract (documented in AGENTS.md §5).

## 2. Database layout

- DB = directory `<datadir>/<name>/` containing `sequences` (FASTA) + BLAST index files
  (`sequences.n*` / `sequences.p*`) built beside it.
- **Header convention**: `>DB~~~GENE~~~ACCESSION~~~RESISTANCE<space>PRODUCT`
  (`~~~` = IDSEP). Fields may be missing; RESISTANCE is `;`-separated, sorted, spaces→`_`.
  RESISTANCE is baked into the ID at DB-build time; **no sidecar metadata is consulted at
  screening time**.
- `--setupdb`: for every subdir with readable `sequences`, run `makeblastdb`; afterwards require
  `sequences.nin` or `sequences.pin` to exist.
- **mol_type heuristic** (replicate exactly): concatenate non-header lines, delete `[AGTC]`
  (case-insensitive); if remaining length > 50% of total → `prot`, else `nucl`. Then
  `makeblastdb -in <path> -title <name> -dbtype <type> -logfile /dev/null`.
- DB introspection: `blastdbcmd -info -db <prefix>`; DBTYPE = `prot` iff output contains
  `total residues` else `nucl`. DBTYPE selects **blastn vs blastx** at screening time.

## 3. Screening pipeline

Per input file (upstream wraps in `bash -c 'set -euo pipefail; ...'`; gapit uses argv lists, no
shell):

```
any2fasta -q -u <file>  |  blastn -task blastn -dust no -perc_identity <minid> \
  -db <datadir>/<db>/sequences \
  -outfmt "6 qseqid qstart qend qlen sseqid sstart send slen sstrand evalue length pident gaps gapopen stitle" \
  -num_threads <threads> -evalue 1E-20 -culling_limit 1 -max_target_seqs 10000
```

- **outfmt 6 fields, exact order (15)**: `qseqid qstart qend qlen sseqid sstart send slen sstrand
  evalue length pident gaps gapopen stitle`. A row with ≠15 columns is a hard error.
- Protein DBs use `blastx -task blastx-fast -seg no` **without `-perc_identity`**, and `--minid`
  is then silently ignored (upstream quirk — reproduce, with a stderr note).
- `any2fasta -q -u` normalizes `.fa/.faa/.gbk/.embl`, gz, bz2 → FASTA on stdout.

## 4. Hit processing (the core algorithm)

For each BLAST row, in order:

1. **Minus-strand normalize**: if `sstrand == "minus"`, swap `sstart`/`send` (subject coords only;
   `qstart`/`qend` are always ascending and never swapped).
2. **Dedup (the only "merge")**: drop the row if `qseqid~qstart~qend` was already seen for this
   input file. First row wins; BLAST emits best hits first. The key **ignores strand**.
   There is **no interval-overlap merging, no gap tolerance, no best-gene choice** — two genes
   overlapping at different query spans are both reported (upstream README caveat). gapit MUST NOT
   add merging to the default path; any future merge mode goes behind a flag, off by default.
3. **Coverage filter**: `pct_cov = 100 * (length - gaps) / slen` (ungapped aligned columns over
   full subject/gene length). Keep iff `pct_cov >= mincov` — comparison on the **unrounded**
   float; display is `%.2f`. (A 79.996% hit displays `80.00` but is discarded. Replicate exactly.)
4. **Identity**: `%IDENTITY = pident` as printed by BLAST (`%.2f`). Never post-filtered.
5. **Subject ID parse**: split `sseqid` on `~~~` → `(database, gene, accession, resistance)`.
   No `~~~` at all → `gene = sseqid`, `accession = ""`, `database = --db` value. Partial fields →
   empty strings.
6. **Product cleanup**: `stitle or "n/a"`; strip all `,` and `\t`; if it contains `~~~`, drop the
   leading whitespace-delimited token (makeblastdb prefixes stitle with the ID; upstream issue #95).

### COVERAGE_MAP (minimap) — replicate this arithmetic exactly

```
WIDTH = 15 - (1 if gapopen > 0 else 0)     # broken maps: 14 boxes + '/'
scale = slen / WIDTH                        # float division
x = int(sstart / scale); y = int(send / scale)
for i in 0 .. WIDTH-1:
    char = '=' if x <= i <= y else '.'
    if gapopen > 0 and i == int(WIDTH/2): append '/' after char
```

Result is always 15 chars. Note the `int()` truncation quirks (e.g. a full-length hit on a long
gene may leave box 0 as `.` since coords are 1-based) — do not "fix" them.

### Row assembly & ordering

| Column | Value |
|---|---|
| FILE | path as given (basename if `--nopath`) |
| SEQUENCE, START, END | `qseqid`, `qstart`, `qend` (query coords) |
| STRAND | `-` if `sstrand == "minus"` else `+` (blastx: always `+`) |
| GENE / DATABASE / ACCESSION / RESISTANCE | `~~~` fields 2/1/3/4 (fallbacks per step 5) |
| COVERAGE | `sstart-send/slen` (subject coords, ascending) |
| COVERAGE_MAP | minimap above |
| GAPS | `gapopen/gaps` |
| %COVERAGE / %IDENTITY | `%.2f` |
| PRODUCT | cleaned stitle |

Rows are sorted by SEQUENCE (lexicographic) then START (numeric) with a **stable** sort, and
emitted per input file after that file finishes. Files are processed sequentially in argument
order; a single header row precedes all output (even if a later file errors — upstream has no
atomicity; gapit buffers per file but preserves row order).

## 5. TSV/CSV output

- Header (unless `--noheader`), printed once:
  `#FILE SEQUENCE START END STRAND GENE COVERAGE COVERAGE_MAP GAPS %COVERAGE %IDENTITY DATABASE ACCESSION PRODUCT RESISTANCE`
- Separator: tab, or `,` with `--csv` (gapit: `--format csv`). Line = `join(sep, fields) + "\n"`.
- stdout = data only; all chatter (Processing/Found N genes/Tips) to stderr.

## 6. Summary mode (`gapit summary`)

Input: ≥1 abricate-format report files. Upstream surface is `abricate --summary`;
gapit exposes it as the `summary` subcommand (`--identity`, `--nopath`, `--format
tsv|csv|json|md`, `--quiet`).

- **Dutch mode**: with exactly 1 input file, matrix rows are keyed by that report's FILE column;
  with >1 files, rows are keyed by input filename (basename applied to labels if `--nopath`).
- First encountered row anywhere is treated as the header map (name→index, later duplicate
  names win — Perl `zip` semantics); lines whose col0 starts with `#` are skipped afterwards.
  A first row WITHOUT `#` (e.g. a noheader report) is header map **and** a data row — quirk
  verified against 1.4.0.
- Duplicate input filenames (compared as given, pre-basename): stderr `WARNING: Skipping
  duplicate file: <name>` + skip. Zero-hit files still appear (non-dutch) with `NUM_FOUND 0`.
- Gene universe = union of all GENE values, sorted lexicographically.
- Output: `#FILE  NUM_FOUND  <gene…>`; cell = each hit's %COVERAGE (or %IDENTITY with
  `--identity`) `;`-joined in file order; absent = `.`. NUM_FOUND = count of **distinct genes**.
- **Rows sort by the as-given key, not the display label** — with `--nopath`, labels can appear
  unsorted (e.g. keys `1dir/zeta.tsv`, `2dir/mid.tsv` print as `zeta.tsv`, `mid.tsv`). Verified.
- The `#`-skip tests raw col0 BEFORE `--nopath` basename-ing. A file key starting with `#` is
  skipped entirely.
- **`[gapit-extension]` divergences** (typed where upstream is silent-undef):
  - separator auto-detected per input file (tab if the first line has one, else comma, else
    tab) — upstream `--csv` must match the reports' format and switches BOTH input and output
    separator; gapit's input parsing is format-agnostic and `--format` controls only the output.
  - no args → UsageError exit 2 (upstream: exit 1); missing/unreadable file → InputError
    `INPUT_NOT_FOUND` exit 5 (upstream: exit 1); non-UTF-8 → InputError `SUMMARY_MALFORMED`.
  - a data row shorter than the mapped GENE/metric columns, or a header map lacking them, →
    InputError `SUMMARY_MALFORMED` (exit 5) with file+line context. Upstream silently treats
    missing cells as empty strings.
  - empty (0-byte) and header-only files are valid: non-dutch they appear with `NUM_FOUND 0`;
    a single one in dutch mode yields the header line only (verified upstream).
- **`gapit.summary/1`** (JSON): `schema`, `tool`, `created_at`, `params` (`metric`:
  `%COVERAGE`|`%IDENTITY`, `nopath`), `genes` (sorted union), `rows` (`file` label,
  `num_found`, `cells`: gene → list of original value strings; absent genes omitted).
  Cells keep the original report strings verbatim. Markdown mirrors this with YAML
  frontmatter and pipe-escaped cells. Introspect via `gapit schema summary`.
- Opt-in parity: `pixi run -e parity summary-parity` byte-diffs gapit vs abricate over the
  committed synthetic fixtures in `tests/data/summary` (dutch, multi, identity, duplicate,
  csv, nopath-key-sort).

## 7. Edge cases & quirks (parity-critical)

- **Input types**: fa/gz/bz2/gbk/embl via any2fasta. Invalid input → pipeline failure → nonzero
  exit (gapit: exit 5 + JSON error envelope).
- **Empty FASTA** → zero hits; header still printed; success exit.
- **Circular contigs**: no special handling (linear).
- **Partial `~~~` headers**: missing trailing fields → empty strings.
- **`--csv` + summary**: summary must be told the separator; mixing breaks upstream — gapit
  detects format per file instead `[gapit-extension]` (RESOLVED 2026-09-17, see §6).
- **Determinism**: stable sort; fixed BLAST params; no timestamps in data payloads. (Upstream
  MOTD/`srand` is stderr-only and dropped in gapit.)
- **blastx**: no sstrand → STRAND `+`; minid unenforced (quirk kept, documented).
- Known upstream caveats we inherit: no mutational resistance; gap reporting incomplete;
  overlapping genes both reported; possible coverage-calculation issues.

## 8. Bundled databases

12 DBs ship in abricate's `db/` (all `nucl` as distributed): `ncbi` (default, AMRFinderPlus core
AMR), `card` (protein-homolog models only; ACC carries `start-end` coords; multi-class
RESISTANCE), `resfinder` (plain ACC; gene has `_copy` suffix), `argannot`, `plasmidfinder`,
`megares` (SNP-confirmation entries excluded), `ecoh`, `vfdb`, `ecoli_vf`, `bacmet2`, `victors`,
`upec_expec_vf`; plus `db/abricate/` — a cd-hit-est recipe, not a DB. Disabled getters:
`ncbibetalactamase`, `serotypefinder`.

gapit v1.0 reads any abricate-format datadir (including abricate's own). A `gapit db fetch`
reimplementation of `abricate-get_db` is post-1.0 (see PLAN.md).

## 9. License note

abricate is GPL-2.0. gapit is a behavioral reimplementation (no Perl code copied); to keep DB
handling and redistribution unambiguous, gapit is licensed GPL-2.0-compatible. Bundled DB content
retains its original upstream licenses.

## 10. Read screening (FASTQ) — gapit extension

abricate cannot screen raw reads; gapit can. Design decisions (2026-09-15):

- **Backend: minimap2 only.** `-x sr` (short/paired reads), `-x map-ont`, `-x map-hifi`.
  bwa/bowtie2 are deliberately excluded: srst2 needs them for SNP-level allele calling, which is
  outside gapit's mission (presence + confidence). One backend covers all read types, and PAF
  output needs no samtools.
- **Invocation**: `gapit screen --r1 R1[,R1b…] [--r2 R2[,R2b…]] --read-type
  sr|map-ont|map-hifi` (comma-separated file lists, one entry per lane; `--r2` count must equal
  `--r1` count — lane i pairs r1[i]/r2[i]) — mutually exclusive with positional contig files.
  Pipeline: one `minimap2 -x <preset> -t <threads> <datadir>/<db>/sequences R1 [R2]` run per
  lane → PAF on stdout; all lanes' rows are aggregated as one sample. The db
  `sequences` FASTA is used directly (minimap2 indexes in memory; no `.mmi` persisted in v1).
- **PAF parsing**: 12 required fields (`qname qlen qstart qend strand tname tlen tstart tend
  nmatch alen mapq`) + optional tags. Keep primary alignments only (`tp:A:P`; records lacking a
  `tp` tag are kept). A row with <12 fields is a hard error.
- **Per-gene aggregation** over a per-base coverage array of length `tlen`:
  `breadth_pct = 100 * covered_bases / tlen` (covered = ≥1 aligned base, union of all primary
  alignments' `tstart..tend`), `mean_depth = sum(per-base depth) / tlen`,
  `reads_mapped` = distinct query names with ≥1 primary alignment on the gene.
- **Presence call**: `present = breadth_pct >= min_breadth`, `--min-breadth` default **90.0**
  (srst2-style). `--minid`/`--mincov` do NOT apply to reads mode.
- **Output**: `--format json` (default in reads mode) emits `gapit.reads/1`; `--format md` the
  Markdown form; `--format tsv|csv` in reads mode is a usage error (exit 2). Genes with zero
  mapped reads are omitted; entries sorted by `breadth_pct` descending, then gene name.
  `gapit.report/1` (contig mode) is unchanged and frozen. `gapit schema reads` introspects the
  new document.
- **Known limits** (document, do not fix in v1): homologous gene families share multi-mapping
  reads — primary-only assignment may misassign closely related alleles; no SNP-level allele
  calling; no per-read identity/MAPQ filtering. Empirically, `minimap2 -x sr` soft-clips ~5 nt
  at each alignment end, so genes under ~100 nt can cap below the 90% breadth threshold — reads
  mode targets normal-length genes (hundreds of nt+); calibrate `--min-breadth` for tiny DBs.
