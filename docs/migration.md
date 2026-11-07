# Migrating from gspec (pre-0.8)

`gspec` is the same software as `gapit` 0.7.x under a new name. The rename
exists because "gspec" collides in literature and code search with
**GAPIT**, the genome-wide association package — every search for the tool
returned the GWAS package instead.

Nothing about the database or output contracts changed. This page lists
exactly what moved, what was added as an alias, and what stayed frozen.

## What changed

| Surface | Before (≤0.7.3) | From 0.8.1 |
|---|---|---|
| PyPI distribution | `gapit` | `gspec` (`gapit` becomes a shim that depends on `gspec`) |
| Command | `gapit` | `gspec` — the `gapit` command still ships and works |
| Version banner | `gspec 0.7.3` | `gspec 0.8.1` |
| `--version --json` `name` | `gapit` | `gspec` (the `gapit.version/1` schema id is unchanged) |
| Datadir env var | `GAPIT_DATADIR` | `GSPEC_DATADIR`; the old name is still read |
| Engine env var | `GAPIT_ENGINE` | `GSPEC_ENGINE`; the old name is still read |
| Home datadir | `~/.local/share/gapit/db` | `~/.local/share/gspec/db`, falling back to the old path when it is the only one present |
| Docs, README, CITATION | `gapit` | `gspec` |

## What deliberately did not change

These are frozen contracts; renaming them would break existing databases,
archived results and pipelines, so they keep the pre-0.8 spelling:

- **JSON schema ids** — `gapit.report/1`, `gapit.typing/1|2|3`,
  `gapit.typing_result/1`, `gapit.cluster/1|2`, `gapit.error/1`, and every
  other `gapit.<document>/<version>` name.
- **Database manifest filename** — `gapit-manifest.json`.
- **Native FASTA header prefix** — `gapit|db=…|gene=…|acc=…|func=…`, and the
  manifest `header_format` value `gapit/v1`.
- **Python import package** — the implementation imports as `gspec`. The
  pre-0.8 `import gapit` keeps working through the `gapit` compatibility
  package (`pip install gapit`), which aliases the same modules, so
  `except gapit.errors.GapitError` still catches what `gspec` raises.
- **Existing database directories** — unchanged on disk and still readable.

## Upgrading

```bash
pip install --upgrade gspec        # or: pip install --upgrade gapit
gspec db path                      # confirm which datadir you are on
gspec db list                      # existing databases still resolve
```

Both commands resolve the same datadir, so a pipeline that calls `gapit` keeps
working while you migrate scripts to `gspec` at your own pace.

To adopt the new home directory explicitly:

```bash
mv ~/.local/share/gapit ~/.local/share/gspec
```

## Timeline

- **0.8.x** — both names ship; old env vars honoured; schemas frozen.
- **1.0** — planned: `gspec.*` schema ids published alongside the `gapit.*`
  ones (aliases, not replacements), and the `gapit` command/import kept as a
  compatibility layer. Nothing is removed without a major version and a
  deprecation cycle.

## Citing

Cite `gspec` (see `CITATION.cff`). For work that used 0.7.x, the old name
still identifies the same software; the Zenodo record for 0.7.3 remains
citable and the 0.8.1 record declares it as a previous version.
