# gapit (transition package)

**`gapit` was renamed to [`gspec`](https://pypi.org/project/gspec/).**

The rename exists because "gapit" collides in literature and code search with
GAPIT, the genome-wide association package, so users looking for this tool
found the GWAS package instead.

## What this package does

Nothing beyond depending on `gspec`. `pip install gapit` now installs the
`gspec` distribution, which provides:

- the **`gspec`** command,
- the **`gapit`** command (kept as an alias, so existing pipelines keep
  working),
- `gspec-mcp` and `gapit-mcp`,
- the **`gapit`** Python import package (`import gapit`) and every frozen
  `gapit.*` JSON schema id.

This distribution ships no console scripts of its own: two distributions
declaring the same script would race on the file in `bin/`.

## What to do

```bash
pip install --upgrade gspec       # or: pip install --upgrade gapit
gspec --version                   # gspec 0.8.0
```

Existing environment variables and data keep working: `GAPIT_DATADIR` is still
read when `GSPEC_DATADIR` is unset, and an existing
`~/.local/share/gapit/db` is still used as the home default. The migration
guide has the full before/after table:
<https://indexofire.github.io/gspec/migration/>

## Version mapping

| `gapit` | Meaning |
|---|---|
| ≤ 0.7.3 | the tool, before the rename |
| 0.8.0 | this transition release: metadata only, installs `gspec>=0.8,<0.9` |

Cite `gspec` (0.8.0 and later). The 0.7.3 Zenodo record remains citable and
identifies the same software.

## Uninstalling

`pip uninstall gapit` removes only this stub; the `gspec` distribution stays
installed (pip does not cascade-dependents), so the `gapit` and `gspec`
commands keep working.
