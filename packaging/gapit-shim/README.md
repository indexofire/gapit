# gapit (transition package)

**`gapit` was renamed to [`gspec`](https://pypi.org/project/gspec/).**

The rename exists because "gapit" collides in literature and code search with
GAPIT, the genome-wide association package, so users looking for this tool
found the GWAS package instead.

## What this package does

It depends on `gspec` and restores the pre-0.9 surface:

- the **`gapit`** and **`gapit-mcp`** commands (they run `gspec`),
- the **`gapit`** Python import package, aliased to the `gspec` tree so
  `except gapit.errors.GapitError` still catches what `gspec` raises.

`gspec` itself provides the `gspec` / `gspec-mcp` commands and the implementation.
The legacy scripts live here, not in `gspec`, because on an upgrade pip installs
`gspec` first and then uninstalls the pre-0.9 `gapit` — which deletes
`bin/gapit` — and this package is installed last, so its files survive.

## What to do

```bash
pip install --upgrade gspec       # or: pip install --upgrade gapit
gspec --version                   # gspec 0.8.0
```

Existing environment variables and data keep working: `GAPIT_DATADIR` is still
read when `GSPEC_DATADIR` is unset, an existing `~/.local/share/gapit/db` is
still used as the home default, and databases written by 0.8
(`gapit-manifest.json`, `gapit|` headers) are read as-is. The migration
guide has the full before/after table:
<https://indexofire.github.io/gspec/migration/>

## Version mapping

| `gapit` | Meaning |
|---|---|
| ≤ 0.7.3 | the tool, before the rename |
| 0.9.0 | this compatibility release: installs `gspec>=0.8.1,<0.10` and keeps `import gapit` + the `gapit` commands |

Cite `gspec` (0.8.0 and later). The 0.7.3 Zenodo record remains citable and
identifies the same software.

## Uninstalling

`pip uninstall gapit` removes only this stub; the `gspec` distribution stays
installed (pip does not cascade-dependents), so the `gapit` and `gspec`
commands keep working.
