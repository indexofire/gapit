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

The pixi `default` and `py312` environments use Python 3.12 and include `gseek>=0.9.1`. Without gseek, `GAPIT_ENGINE=gseek` raises a
typed `MISSING_DEPENDENCY` error; `legacy` and all existing commands are
unaffected.

## Why legacy is retained

This is a deliberate dual-path design:

- `legacy` is the default for every installation;
- every migrated path keeps its original implementation and remains usable;
- `shadow` compares both paths and logs coarse differences while returning the
  legacy result;
- legacy is the parity baseline and the fallback when gseek is unavailable;
- database building, typing, formats, MCP, and CLI behavior are **not** routed
  through gseek.

## Verification

The bridge has unit coverage for `GAPIT_ENGINE`, coordinate/strand mapping,
protein task selection, PAF/NM/cs tags, dispatch, and shadow fallback. Real
parity tests build tiny nucleotide/protein BLAST databases and a minimap2
target, then compare legacy and gseek rows exactly; they skip when gseek or the
external binaries are missing.
