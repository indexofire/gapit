# gseek engine (optional)

gapit can route its alignment subprocesses through
[gseek](https://github.com/indexofire/gseek), but the existing subprocess
implementation remains the default and is always available.

## Switch

```bash
GAPIT_ENGINE=legacy   # default: existing blastn/blastx/minimap2 calls
GAPIT_ENGINE=gseek    # gseek plans and parses the same alignments
GAPIT_ENGINE=shadow   # run both, log coarse differences, return legacy
```

`GAPIT_ENGINE` is read per call. An unknown value fails before any alignment
runs. `shadow` is a verification mode, not a third production path: it always
returns the legacy result.

## What is routed

| Path | legacy | gseek |
| --- | --- | --- |
| Contig screening BLAST | `blastn`/`blastx` + tabular parser | `gseek.align(backend="blast")` |
| Protein screening | `blastx -task blastx-fast` | `gseek.align(backend="blast", task="blastx-fast")` |
| Reads screening | streaming `minimap2` PAF | `gseek.align(backend="minimap2")` |
| Cluster screening | `minimap2 -x asm20 --cs` | `gseek.align(backend="minimap2", preset="asm20", cs=True)` |

gseek mode keeps gapit's thresholds, coordinate conventions, and output
contracts; it only changes who plans the command and parses the rows. In gseek
mode minimap2 PAF rows are materialized per lane by gseek's API, whereas the
legacy path streams them. gseek returns a deterministic sorted record view, so
raw PAF emission order can differ from minimap2; gapit's reads/cluster
consumers aggregate and sort, so the final reports are unchanged. The bridge
also preserves two consumer-specific compatibility cases: VFDB's duplicate
subject IDs are kept by gseek, and CARD's >255-character target names use a
direct FASTA target (`index="none"`) because minimap2's `.mmi` index replaces
such names with numeric sequence indices.

## Install

gapit requires Python >= 3.12, and every supported environment includes
gseek:

```bash
pip install "gapit[gseek]"
```

The pixi `default` environment (with the `py312` feature) uses Python 3.12 and includes `gseek>=0.10.2`. Without gseek, `GAPIT_ENGINE=gseek` raises a
typed `MISSING_DEPENDENCY` error; `legacy` and all existing commands are
unaffected.

## Why legacy is retained

The switch itself uses `gseek.compat.Engine`, so the legacy/gseek/shadow
selection and fallback are shared with the other consumers; gapit keeps only
its own legacy implementation, gseek call, and result semantics.

This is a deliberate dual-path design:

- `legacy` is the default for every installation;
- every migrated path keeps its original implementation and remains usable;
- `shadow` compares both paths and logs coarse differences while returning the
  legacy result;
- legacy is the parity baseline and the fallback when gseek is unavailable;
- database building, typing, formats, MCP, and CLI behavior are **not** routed
  through gseek.

## Performance characteristics

The two engines optimize different parts of the pipeline, so the faster one depends on the
database kind:

- **Gene BLAST (contig mode).** Legacy invokes `blastn`/`blastx` against the database's
  existing BLAST index, once per input file. gseek maintains its own content-addressed
  BLAST index under `~/.cache/gseek` and still plans/parses each query. gseek >= 0.10.2
  probes only the tools the task needs; the bridge additionally reuses one
  `ToolResolver`/`Executor`/`IndexStore` per process (a `PATH` or `GSEEK_*` change
  starts a fresh session), so neither layer re-probes every BLAST binary per file. When
  the installed gseek exposes the precomputed-ref/external-DB hooks, the bridge hashes the
  target once and points gseek at gapit's existing BLAST DB; older gseek versions keep
  their own index. If gseek also exposes `align_many` and every input is FASTA, gapit
  batches them into one invocation (split into at most `--jobs` chunks), which amortizes
  process/index setup across queries. On a gseek build with the newer batch options,
  gapit skips query-input hashing and gseek's canonical sort because it folds and orders
  hits itself. Otherwise the bridge stays per-file.
- **Reads and cluster minimap2.** The gseek bridge uses `iter_align` when available, so
  PAF records flow straight into gapit's breadth/cluster math instead of being
  materialized by `gseek.align`.
- **Cluster minimap2 (assembly mode).** A cluster build (and `db setup`) writes
  `<db>/sequences.asm20.mmi`; legacy uses that prebuilt index when it is fresh, so it does
  not reindex the locus FASTA per file. gseek keeps its own content-addressed `.mmi` under
  `~/.cache/gseek`, so both paths can reuse an index across files; gseek's per-call
  planning still costs something. Without a prebuilt index legacy falls back to the FASTA
  target and rebuilds it in memory for every input. On a database with a fresh index legacy
  can be competitive with, or faster than, gseek; without one, gseek usually wins
  multi-file runs.
- **Cache warmth.** The first gseek run pays the index build; later runs reuse the cache.
  Index keys include target content, semantic parameters, adapter version, and tool
  versions. `GSEEK_CACHE_DIR` relocates the cache; it is safe to delete when diagnosing
  stale indexes.

`--jobs` amplifies these effects: each worker has its own aligner process, while the
gseek session and index cache are shared. Benchmark with warm caches and the same
`--threads`/`--jobs` when comparing engines.

## Verification

The bridge has unit coverage for `GAPIT_ENGINE`, coordinate/strand mapping,
protein task selection, PAF/NM/cs tags, dispatch, and shadow fallback. Real
parity tests build tiny nucleotide/protein BLAST databases and a minimap2
target, then compare legacy and gseek rows exactly; they skip when gseek or the
external binaries are missing.
