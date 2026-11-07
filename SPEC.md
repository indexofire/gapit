# SPEC.md — gspec legacy-table and compatibility contract

> **Scope.** This document is the frozen interoperability contract for the
> abricate-shaped surface gspec still consumes and emits: `~~~` database
> headers, the 15-column screen table, the value rules behind its fields, and
> the table shapes `gspec typing` / `gspec summary` accept. It is **not** the
> abricate CLI specification and it is not the only gspec contract.
>
> Current gspec features — subcommands, reads, cluster screening, typing,
> database providers/bundles, gseek, JSON documents — are specified by the
> docs and schemas linked in §7. The abricate parity harness remains a
> regression reference for this legacy surface; from v0.5.0 it is not a
> design goal.
>
> **Frozen** means field names/order/units and the named arithmetic below do
> not change without the compatibility decision in AGENTS.md §5/§10. New
> gspec surfaces may diverge from abricate by design.

## 1. Scope and compatibility policy

| Surface | Status |
|---|---|
| `~~~` database headers (§2) | frozen read contract |
| Legacy contig BLAST path and its field values (§3–§4) | frozen |
| 15-column screen table (§5) | frozen |
| `gspec typing` / `gspec summary` legacy-table inputs (§6) | frozen shape |
| Current CLI, subcommands, exit codes, schemas | gspec-owned; AGENTS.md §5 |
| Reads, clusters, typing rules, databases/providers, gseek | gspec-owned; §7 links |

Where this file and a linked current-feature doc disagree about a gspec-only
feature, the linked doc wins. For §2–§6, this file wins for the frozen
interop surface.

The parity harness (`pixi run -e parity parity` / `summary-parity`) checks
the legacy surface against real abricate. It is opt-in, regression-only, and
does not constrain new gspec features.

## 2. Legacy database layout and `~~~` headers `[frozen]`

- A legacy database is a directory `<datadir>/<name>/` containing a
  nucleotide or protein FASTA named `sequences`, with the BLAST index
  (`sequences.n*` / `sequences.p*`) built beside it.
- **Header convention**:
  `>DB~~~GENE~~~ACCESSION~~~RESISTANCE<space>PRODUCT` (`~~~` = IDSEP).
  Fields may be missing. `RESISTANCE` is `;`-separated, sorted, and spaces
  become `_`. For legacy databases the RESISTANCE value is baked into the ID;
  no sidecar metadata is consulted at screening time.
- Decode fallbacks: split on `~~~`; no separator at all means
  `gene = sseqid`, `accession = ""`, `database = --db`; a partial split fills
  missing trailing fields with empty strings.
- Legacy indexing runs `makeblastdb` with an explicit `-dbtype`, selected by
  the abricate mol_type heuristic: concatenate non-header sequence, delete
  `[AGTC]` case-insensitively; if the remainder is >50% of the total it is
  `prot`, else `nucl`. Current `gspec db setup` / `gspec db build` and the
  native `gapit/v1` format are documented in `docs/databases.md`; gspec reads
  both native and legacy datadirs.
- DBTYPE introspection uses `blastdbcmd -info`; `prot` iff the output
  contains `total residues`, else `nucl`. DBTYPE selects `blastn` vs `blastx`.

## 3. Legacy contig screening path `[frozen]`

These are the defaults and the BLAST invocation behind the legacy table.
Reads and cluster modes are separate gspec contracts (§7).

- Defaults for the legacy contig path: `--minid 80`, `--mincov 80`,
  `--threads 1`. `--minid` is enforced only through `blastn -perc_identity`;
  `--mincov` is a post-filter on the unrounded coverage float.
- BLAST invocation:

  ```
  blastn -task blastn -dust no -perc_identity <minid> \
    -db <datadir>/<db>/sequences \
    -outfmt "6 qseqid qstart qend qlen sseqid sstart send slen sstrand evalue length pident gaps gapopen stitle" \
    -num_threads <threads> -evalue 1E-20 -culling_limit 1 -max_target_seqs 10000
  ```

  The outfmt field order is exact and has 15 fields; a row with a different
  count is a hard error.
- Protein DBs use `blastx -task blastx-fast -seg no` **without**
  `-perc_identity`; the legacy path silently ignores `--minid` and emits a
  stderr note. This quirk is retained for compatibility.
- Input normalization is native (`seqconvert.py`, any2fasta `-q -u`
  semantics): fa/fq/gbk/embl, gz/bz2 → FASTA. gspec re-wraps at 60 columns and
  reads universal newlines; these are blast-invisible. Current screening
  options (`--fofn`, `--jobs`, CSV/JSON/MD, progress) are in
  `docs/screen.md`.
- The legacy path checks that `blastn` is ≥ 2.2.30 before running; current
  environments require modern BLAST+ anyway.

## 4. Hit processing and COVERAGE_MAP `[frozen]`

For each BLAST row, in order:

1. **Minus-strand normalize**: if `sstrand == "minus"`, swap `sstart`/`send`
   (subject coords only; query coords are already ascending).
2. **Dedup (the only "merge")**: drop a row if `qseqid~qstart~qend` was already
   seen for this input file; first row wins; the key ignores strand. There is
   no interval-overlap merging and no best-gene choice in the default path.
   `--merge-fragments` is an explicit gspec extension
   (`docs/screen.md`), not part of this frozen path.
3. **Coverage filter**:
   `pct_cov = 100 * (length - gaps) / slen` (ungapped aligned columns over
   full subject length). Keep iff `pct_cov >= mincov`, compared on the
   unrounded float; display is `%.2f`. A 79.996% hit displays `80.00` but is
   discarded.
4. **Identity**: `%IDENTITY = pident` as printed by BLAST (`%.2f`), never
   post-filtered.
5. **Subject ID parse**: legacy `~~~` rules above; native `gapit/v1` headers
   decode through `docs/databases.md` and are not part of the legacy grammar.
6. **Product cleanup**: `stitle or "n/a"`; strip `,` and `\t`; if the product
   contains `~~~` (or the subject is a gspec-tagged seqid), drop the leading
   whitespace-delimited token (makeblastdb prefixes stitle with the ID).

### COVERAGE_MAP (minimap)

```
WIDTH = 15 - (1 if gapopen > 0 else 0)     # broken maps: 14 boxes + '/'
scale = slen / WIDTH                        # float division
x = int(sstart / scale); y = int(send / scale)
for i in 0 .. WIDTH-1:
    char = '=' if x <= i <= y else '.'
    if gapopen > 0 and i == int(WIDTH/2): append '/' after char
```

The result is always 15 chars. The `int()` truncation quirks are intentional
and must not be "fixed".

### Row assembly and ordering

| Column | Value |
|---|---|
| FILE | path as given (basename with `--nopath`) |
| SEQUENCE, START, END | `qseqid`, `qstart`, `qend` (query coords) |
| STRAND | `-` if `sstrand == "minus"`, else `+` (blastx: always `+`) |
| GENE / DATABASE / ACCESSION / RESISTANCE | `~~~` fields 2/1/3/4 (fallbacks above) |
| COVERAGE | `sstart-send/slen` (subject coords, ascending) |
| COVERAGE_MAP | the minimap above |
| GAPS | `gapopen/gaps` |
| %COVERAGE / %IDENTITY | `%.2f` |
| PRODUCT | cleaned stitle |

Rows sort by SEQUENCE lexicographically, then START numerically, with a
stable sort. Files are processed in argument order; `gspec screen --jobs N`
may screen concurrently but stdout stays in input order. The TSV header is
printed once, and completed files stream out as they finish; a later failure
leaves the already-emitted prefix plus the typed error envelope.

## 5. 15-column table `[frozen]`

- Header unless `--noheader`:
  `#FILE SEQUENCE START END STRAND GENE COVERAGE COVERAGE_MAP GAPS %COVERAGE %IDENTITY DATABASE ACCESSION PRODUCT RESISTANCE`
- Separator: tab, or comma via `--format csv`. The legacy abricate `--csv`
  flag is not a gspec CLI option; `--format csv` is the gspec spelling.
- One line per hit: `join(separator, fields) + "\n"`. stdout carries data
  only; diagnostics, progress and the `Found N genes` note go to stderr.
- `--nopath` applies to the FILE column. `RESISTANCE` is the display/resistance
  field; for native databases it carries the locked function vocabulary
  documented in `docs/databases.md`.
- `#`-prefixed lines after the header are skipped by table consumers; this
  allows concatenated screen tables.
- JSON and Markdown are gspec extensions, not this 15-column contract; their
  schemas and fields are in `docs/outputs.md` and `gspec schema report`.

## 6. Consuming legacy tables `[frozen shape]`

### `gspec typing`

`gspec typing` reads TSV or CSV screen tables from files or stdin. The contig
shape needs `FILE`/`#FILE`, `GENE`, `DATABASE`, `%IDENTITY`, and `%COVERAGE`;
the reads shape needs `#SAMPLE`, `GENE`, `DATABASE`, and `BREADTH%`. A table
may contain rows for several FILE values but only one DATABASE value; mixed
DATABASE values are a typed error. `#`-prefixed continuation header rows are
skipped. Typing uses the values as given; it does not re-run BLAST. Current
rule semantics, reads input details, and the output document are in
`docs/typing.md` and `gspec schema typing_result`.

### `gspec summary`

`gspec summary` accepts one or more abricate-format report tables (TSV or CSV,
auto-detected per file). The first line is the header map; later lines whose
first cell starts with `#` are skipped. Rows are keyed by FILE (dutch mode) or
input filename (multi-file); duplicate input filenames are warned and skipped.
The gene universe is the sorted union of GENE values; `NUM_FOUND` counts
distinct genes; cells are `;`-joined original metric strings in file order
(`.` when absent). Current modes, formats, JSON, and Markdown are in
`docs/summary.md`; the JSON document is `gspec schema summary`.

## 7. Current gspec contract map

| Surface | Current contract |
|---|---|
| CLI, exit codes, error envelope | `AGENTS.md` §5, `gspec --help`, `gspec schema error` |
| Contig screening, formats, fragments | `docs/screen.md`, `docs/outputs.md` |
| Reads (FASTQ/FASTA, minimap2) | `docs/reads.md`, `gspec schema reads` / `reads2` |
| Cluster databases (`kind: cluster`) | `docs/screen.md#cluster-databases-kind-cluster` |
| Typing (`gspec typing`) | `docs/typing.md`, `gspec schema typing` / `typing_result` |
| Databases, providers, bundles, headers | `docs/databases.md`, `docs/datadir.md` |
| gseek engine | `docs/gseek-engine.md` |
| MCP | `docs/mcp.md` |
| Output schemas | `gspec schema <name>` (`gapit.schema` registry) |

## 8. Change policy

- A change to §2–§6 (legacy header, field value, table shape, typing/summary
  input) is a compatibility change: update this file, the relevant
  tests/goldens, and the parity harness.
- A change to a current gspec-only feature (new subcommand, reads/cluster
  behavior, typing rule, provider, bundle, gseek) updates its linked doc and
  tests; it does **not** require a SPEC.md rewrite.
- JSON field or schema changes follow AGENTS.md §5/§10: bump the schema
  version and record the change; never rename or retype a field in a minor
  bump.
- The parity harness is a regression reference for the legacy surface, not a
  release gate for new features.

## 9. License and provenance

gspec is MIT-licensed; abricate is GPL-2.0. gspec is an independent behavioral
reimplementation and copies no Perl code. Legacy database content retains its
upstream license. Bundled, provider, and native-database provenance rules live
in `AGENTS.md` §6 and `docs/databases.md`; this file is documentation, not
code.
