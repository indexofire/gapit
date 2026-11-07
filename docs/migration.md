# Migrating from gapit

`gspec` is the same software as `gapit` under a new name. The rename exists
because "gapit" collides in literature and code search with **GAPIT**, the
genome-wide association package — every search for this tool returned the GWAS
package instead.

The rename happened in two steps:

| Release | Renamed | Left alone |
|---|---|---|
| **0.8** | project name, PyPI distribution, CLI banner, docs, bundled vendor strings | the `gapit` command, `import gapit`, every on-disk name |
| **0.9** | every on-disk contract (schema ids, manifest filename, FASTA header prefix) | the legacy spellings, which are still read |

## What changed

| Surface | Before (≤0.8.1) | From 0.9.0 |
|---|---|---|
| PyPI distribution | `gapit` | `gspec` (the `gapit` compatibility package depends on it) |
| Command | `gapit` | `gspec`; `gapit` still works through the compatibility package |
| Version banner | `gapit 0.8.1` | `gspec 0.9.0` |
| `--version --json` | `gapit.version/1`, `name: gapit` | `gspec.version/1`, `name: gspec` |
| JSON schema ids | `gapit.<doc>/<n>` | `gspec.<doc>/<n>` |
| Manifest filename | `gapit-manifest.json` | `gspec-manifest.json` |
| FASTA header prefix | `gapit\|db=…` | `gspec\|db=…` |
| Manifest `header_format` | `gapit/v1` | `gspec/v1` |
| Datadir env var | `GAPIT_DATADIR` | `GSPEC_DATADIR` (old name still read) |
| Engine env var | `GAPIT_ENGINE` | `GSPEC_ENGINE` (old name still read) |
| Home datadir | `~/.local/share/gapit/db` | `~/.local/share/gspec/db`, falling back to the old path when it is the only one present |
| Python import | `import gapit` | `import gspec`; `import gapit` via the compatibility package |

## What still works unchanged

Every pre-0.9 spelling is **accepted on read**, and every writer emits the
current one:

- **Schema ids** — a document whose `"schema"` is `gapit.typing/2`,
  `gapit.report/1`, `gapit.cluster/2`, … validates and is re-emitted with the
  `gspec.*` id. Archived JSON stays readable; nothing needs rewriting.
- **Manifest filename** — `gapit-manifest.json` is found when
  `gspec-manifest.json` is absent. A database installed under 0.8 keeps
  working as-is.
- **FASTA header prefix** — `gapit|db=…|gene=…|acc=…|func=…` still decodes;
  new builds write `gspec|…`.
- **`header_format`** — `gapit/v1` reads as `gspec/v1`.
- **Environment** — `GAPIT_DATADIR` / `GAPIT_ENGINE` are read when the
  `GSPEC_*` names are unset, and a pre-0.9 home datadir is used when the new
  one does not exist.

Nothing forces a migration. To normalise an existing datadir, rebuild the
databases (`gspec db setup --all` reindexes in place; `gspec db build`
rewrites the manifest and headers) or move the directory:

```bash
mv ~/.local/share/gapit ~/.local/share/gspec
```

## Upgrading

```bash
pip install --upgrade gspec        # new name only
pip install --upgrade gapit        # keeps import gapit + the gapit command
gspec db path                      # confirm which datadir you are on
gspec db list                      # existing databases still resolve
```

Both commands read the same datadir and the same databases, so a pipeline that
calls `gapit` keeps working while scripts migrate to `gspec` at their own pace.

## Timeline

- **0.8** — brand rename; on-disk names unchanged.
- **0.9** — on-disk rename; legacy spellings still read.
- **1.0** — planned: the legacy read aliases stay; nothing is removed without
  a major version and a deprecation cycle.

## Citing

Cite `gspec` (see `CITATION.cff`). Work done under 0.x can cite the release it
used; the Zenodo records identify the same software.
