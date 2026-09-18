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
gapit-native DB construction and acquisition (`gapit db`, `gapit/v1` headers) is specified
in §11.

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

## 11. gapit-native databases — [gapit-extension]

gapit builds its own databases rather than only consuming abricate-built datadirs. A native
DB keeps a typed truth source (`records.jsonl`) and generates every screening artifact from
it; sequence headers are tagged + percent-encoded (`gapit/v1`) instead of `~~~`. Design
locked 2026-09-17. **abricate cannot read gapit-built DBs** (accepted trade-off); gapit
reads both formats — no `gapit|` prefix → legacy rules (§4 step 5), untouched.

### Header grammar (`gapit/v1`)

`>gapit|db=<v>|gene=<v>|acc=<v>|func=<v><space>PRODUCT`; fixed key order `db gene acc func`
(writers emit exactly this order; decode keys by name — order-agnostic, duplicates
rejected); empty `acc`/`func` values allowed.

- **Percent-encoding** (values only; lowercase hex; every other char verbatim):

  | raw | encoded | raw | encoded |
  |---|---|---|---|
  | `%` | `%25` | space | `%20` |
  | `\|` | `%7C` | tab | `%09` |
  | `=` | `%3D` | LF | `%0A` |
  | `;` | `%3B` | CR | `%0D` |

- **function**: the list is joined with `;` first, then encoded as one value → decodes to
  the abricate display string `a;b`.
- **encode total / decode strict**: encode never raises. Decode raises `DatabaseError
  HEADER_MALFORMED` (exit 4; context `{seqid, reason}`; machine-stable snake_case reasons —
  `missing_gene`, `missing_db`, `missing_acc`, `duplicate_key:<k>`,
  `invalid_percent_escape:<k>`, `segment_without_key`, `empty_key`) on: a `%` not followed
  by two hex digits (checked in EVERY segment, unknown keys included), missing/empty
  `gene`, a missing `db`/`acc` KEY, duplicate key, segment without `=`, empty key. The
  `db`/`gene`/`acc` keys are required-present (gapit answers presence/absence — identity
  fields must be explicit); an empty `acc` VALUE is allowed (unpublished sequences) and an
  empty `db` value falls back to `default_db` (symmetric with legacy); the `func` key is
  intentionally optional (functional annotation is incidental metadata) — absent `func`
  decodes to `""`. Unknown keys are skipped by name (forward compat). Uppercase-hex
  escapes decode fine (strict-write, lenient-read).
  The same error fires at screening time on a malformed tagged `sseqid`/`tname`.
- **round-trip evidence** (BLAST 2.15+, minimap2): tagged seqids survive makeblastdb →
  blastn `sseqid` and minimap2 `tname` byte-identical; 727-char ids pass untruncated; BLAST
  does not reinterpret pipe ids; `blastdbcmd -info` (§2 DBTYPE introspection) works.
- **product cleanup extension** (issue-#95 quirk translated): stitle is always
  `"<seqid> <product>"`, so §4 step 6's leading-token strip also fires when the seqid is
  gapit-tagged (`~~~` in product **or** tagged seqid); a bare tagged seqid in stitle stays
  verbatim.

### Directory layout

```
<datadir>/<name>/
  records.jsonl         truth source: one frozen Record JSON object per line (LF, UTF-8)
  sequences             GENERATED projection: tagged headers, 60-col wrap, atomic write
  sequences.n*|p*       BLAST index, explicit dbtype (§2 mol_type heuristic not consulted)
  sequences.mmi         minimap2 -d index (nucl only; prot builds skip it silently)
  gapit-manifest.json   gapit.manifest/1 provenance sidecar, written LAST (certifies build)
```

- **Record fields** (frozen): `db`, `gene`, `sequence`, `accession` (default `""`),
  `function` (tuple, default `()`), `product` (default `"n/a"`), `source_id` (default
  `""`).
- **gapit.manifest/1 fields**: `schema`, `name`, `source_urls`, `fetched_at` (ISO-8601 UTC,
  seconds), `sha256` (of `sequences`; 64 hex, stored lowercase), `n_records`, `dbtype`
  (explicit `nucl|prot`), `header_format` (`gapit/v1`), `upstream_version`, `tool`,
  `makeblastdb_version`, `minimap2_version` (the `.mmi` reuse gate, below).
- `records.jsonl` + the manifest are FILE contracts: never on stdout, not in
  `gapit schema`. Rebuilds are deterministic — same records + inputs → byte-identical
  `sequences` and manifest (tool versions come from the environment, `fetched_at` is an
  input; nothing wall-clock).

### Build pipeline

```
fetch → transform → normalize → dedupe → sort → records.jsonl → sequences → self-check
  → sha256 → makeblastdb (explicit dbtype) → .mmi (nucl only) → manifest (written LAST)
```

- **normalize** per record: uppercase; nucl `[^AGCT]`→`N`, prot `[^A-Z]`→`X`; function
  sorted, then whitespace runs → `_` per class (upstream load_fasta semantics).
- **dedupe** on the exact normalized sequence, first wins (dropped count → stderr note);
  duplicate gene NAMES are allowed and kept (upstream keeps them too).
- **sort** by `gene` (lexicographic) → `records.jsonl`.
- **self-check**: every generated header must decode back to its record — in order, counts
  included; product is NOT compared (display field, not truth). Any failure stops the build
  before anything is indexed.
- Upstream `is_full_gene` is a no-op (map result discarded) — not ported.
- **typed errors** (DatabaseError, exit 4): `DB_ALREADY_EXISTS` (manifest present, no
  `--force`), `DOWNLOAD_FAILED`, `PROVIDER_EMPTY`, `PROVIDER_INVALID`, `BUILD_INVALID`
  (product line break / empty sequence / zero records), `BUILD_SELF_CHECK_FAILED`,
  `MMI_BUILD_FAILED`.

### Providers

| name | dbtype | source | transform (distilled) |
|---|---|---|---|
| `ncbi` | nucl | `https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/` (`AMR_CDS.fa` + `ReferenceGeneCatalog.txt`) | 7-field id `pi\|acc\|fp\|fn\|gene\|fam\|prod` (later pipes fold into prod); keep iff `fp==fn==1`; TSV row keyed `refseq_nucleotide_accession` (0-based col 10) required, scope `core`, type/subtype `AMR`; ABX = subclass `/`-split; ACC = TSV acc (`+.1` if unversioned); DESC = prod `_`→space |
| `card` | nucl | `https://card.mcmaster.ca/latest/data` (tar.bz2 → `card.json`) | `protein homolog model` only; `model_param.snp` → PROVIDER_INVALID; dna = first sorted `model_sequences` key; ARO `Drug Class` → ABX (strip ` antibiotic`, spaces→`_`); ID = model_name ws→`_`; ACC = `accession:fmax-fmin` (swapped iff minus strand); DESC = ARO_description |
| `resfinder` | nucl | `https://bitbucket.org/genomicepidemiology/resfinder_db/get/HEAD.zip` | `phenotypes.txt` keyed by gene (`^(.*?)_\w+$`, last row wins); ABX = col2 `,`-split minus `unknown\|notes\|none` (filter before trim); issue-#62 glued-`>` repair on raw text; `.fsa` id `^(.*?)_(\d+)_(\S+)$` → ID = `{id}_{copy}`; ACC; DESC = base gene |
| `argannot` | nucl | `https://www.mediterranee-infection.com/wp-content/uploads/2019/09/ARG-ANNOT_NT_V6_July2019.txt` | strip `\` file-wide first; `:`-split → ID = x0, ACC = `x1:x2`, DESC = `""` |
| `plasmidfinder` | nucl | `https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/HEAD.zip` | `*.fsa`; `^(.*)_(([A-Z]+\|NC_)\d+(\.\d+)?)$` keeps the copy number ON the gene; ACC = accession; DESC = original full id; non-matching id kept verbatim |
| `megares` | nucl | `https://www.meglab.org/downloads/megares_v3.00.zip` | `\|`-split (limit 6); non-empty 6th field (RequiresSNPConfirmation) → skip; ID = group, ACC = MEG id, DESC = `x1:x2:x3:x4`; no ABX |
| `ecoh` | nucl | `https://raw.githubusercontent.com/katholt/srst2/master/data/EcOH.fasta` | `__`-split → gene = x2 (allele field); `;`-split desc → ACC = first, DESC = rest joined with space |
| `vfdb` | nucl | `http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz` | id `^(\w+)\(\w+\|(\w+)(\.\d+)?\)$` → ACC = $2 (version suffix dropped); first `(...)` of DESC renames the gene; DESC otherwise verbatim |
| `ecoli_vf` | nucl | `https://github.com/phac-nml/ecoli_vf/raw/master/data/repaired_ecoli_vfs_shortnames.ffn` | id `^(\w+)(?:\((.*?)\))?$` → ACC = `$2 or $1`; strip trailing `[strain]` groups; paren-desc overrides gene |
| `bacmet2` | **prot** | `http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta` | `\|`-split, ≥4 fields (5th UniProt name ignored); ID = `x1-x0` (reversed); ACC = `x2:x3`; DESC = description or gene |
| `victors` | nucl | `http://phidias.us/victors/downloads/gen_downloads.php` + `gen_downloads_protein.php` | `.faa` headers → gi map (ACC, DESC); `.ffn` id `gi\|N:s-e` → ACC = map or `gi\|N:s-e`; DESC = map desc cut at ` [` or `hypothetical protein`; gene = full `.ffn` id |
| `upec_expec_vf` | nucl | `https://raw.githubusercontent.com/FordeGenomics/ST167_Code/refs/heads/main/UPEC-ExPEC_VF/UPEC_ExPEC_VF.tsv` | TSV keyed `Gene name` (first row wins); ID = gene; ACC = `Accession/Source:Begin-End`; DESC = Description; ABX = [Class]; SEQ = Sequence |

Known divergences vs `abricate-get_db` 1.4.0: records whose ids fail a provider regex are
**skipped** (upstream splices stale/undef perl capture variables into garbage fields);
resfinder + plasmidfinder fetch Bitbucket `HEAD.zip` instead of `git clone` (no git
dependency). Transform semantics are distilled from the perl (read-only) and pinned
per-provider by offline synthetic-fixture unit tests; real-network fetch is the runtime
path and is not exercised by the offline suite.

### CLI (`gapit db`)

- `gapit db fetch NAME [--datadir D] [--force] [--quiet]` — full pipeline into
  `<datadir>/NAME`; stdout = one-line JSON receipt `{db, records, dbtype, destination}`;
  unknown NAME → UsageError exit 2 (registry lookup precedes datadir resolution).
- `gapit db list [--datadir D] [--json]` — `PROVIDER STATUS DBTYPE DESCRIPTION` table
  (`installed (N)` iff `<datadir>/<name>/gapit-manifest.json` exists, else `available`),
  or `gapit.dblist/1` JSON (`providers[]`: `name`, `description`, `dbtype`, `installed`,
  `records` omitted when unset). Module-local schema, deliberately NOT in `gapit schema`.
- `gapit db install SOURCE --sha256 HASH --output TARGET` — verified LOCAL-FILE install
  only (no network, no providers, no archives): streaming SHA256, atomic replace after the
  digest verifies; stdout receipt `{destination, sha256}`; mismatch → InputError
  `CHECKSUM_MISMATCH` exit 5.

### Reads mode: `.mmi` reuse

- `sequences.mmi` is reused iff it exists, a readable manifest sits beside it, and
  `manifest.minimap2_version` equals the installed `minimap2 --version` (version-only gate;
  the `sequences` sha256 is deliberately not re-hashed). Any miss → **silent FASTA
  fallback** (the §10 in-memory path — a performance fallback, not an error); a genuinely
  missing minimap2 surfaces as the typed DependencyError from the real invocation. Legacy
  datadirs carry no manifest → always the FASTA path. This supersedes §10's "no `.mmi`
  persisted in v1" for gapit-native DBs.
- Reads-mode `tname` decodes through the same codec (tagged or `~~~`), so native DBs yield
  proper gene fields in `gapit.reads/1`; prot DBs build no `.mmi` (minimap2 is
  nucleotide-only).

### Function vocabulary (gapit/v1 `func` values)

Each native provider populates `Record.function` from a locked vocabulary (values flow
to the `func=` header key and onward to the outputs):

| provider(s) | `func` value(s) |
|---|---|
| `ncbi`, `resfinder`, `argannot`, `card` | the source's antibiotic classes (the abricate ABX values) |
| `megares` | the MEGARes class field |
| `ecoh` | allele-derived `H-antigen` / `O-antigen` / antigen fallback |
| `vfdb`, `ecoli_vf`, `victors`, `upec_expec_vf` | `virulence` |
| `plasmidfinder` | `replicon` |
| `bacmet2` | `biocide` (BacMet spans biocides + metals; the id carries no class) |

The TSV `RESISTANCE` / JSON `resistance` output names are frozen and carry these
functional categories for native DBs.

### Bundled snapshots

- card and vfdb ship as **bundled snapshots** inside the wheel
  (`src/gapit/data/snapshots/<name>.tar.gz`, Wave G): `gapit db fetch NAME` installs
  them with zero network. All other providers keep the upstream fetch;
  `--from-source` forces the upstream download even when a snapshot exists, and a
  missing/unresolvable archive falls back to the network path silently.
- Archive layout (frozen): `<name>.tar.gz` containing exactly `records.jsonl` +
  `gapit-manifest.json` from an installed db dir. Snapshots carry
  **post-normalize records** — never BLAST/minimap2 indexes (index bytes are
  BLAST-version-sensitive; a local rebuild from records is deterministic and
  fast). Archives are built deterministically: sorted entry names, PAX format,
  gzip mtime 0 → byte-identical rebuilds. The installed manifest rebuilds
  `fetched_at`/`sha256`/tool versions locally; only `upstream_version` is
  inherited from the archived manifest.
- Bare `gapit db fetch` (no NAME) installs the default set `("card", "vfdb")` in
  order, one JSON receipt line per db on stdout. Adding a future bundled DB =
  dropping a `<name>.tar.gz` into the snapshots dir + setting `snapshot=` on the
  provider (one line).
