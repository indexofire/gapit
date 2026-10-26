# Contributing to gapit

gapit is a typed, test-first Python project. Start with `AGENTS.md` for the
repository rules and `.github/workflows/ci.yml` for the exact gates.

## Development setup

Use [pixi](https://pixi.sh) (never pip/conda directly in the dedicated env):

```bash
pixi install
pixi run lint
pixi run fmt-check
pixi run typecheck
pixi run test
```

External BLAST+ and minimap2 binaries arrive through the pixi environment.
The optional `gseek` engine is covered by the default environment; legacy
(`GAPIT_ENGINE=legacy`) is always available.

## Every change

- Keep `lint`, `fmt-check`, `typecheck`, and `test` green.
- Add a regression test for every bug fix; use table-driven tests for edge
  cases (delimiters, malformed tables, empty inputs).
- Preserve deterministic output: stable sorts, LF endings, UTF-8, no wall-clock
  timestamps inside data payloads.
- Keep stdout data-only; diagnostics go to stderr.
- Public contracts (CLI flags, exit codes, JSON schemas) need a schema/version
  note in `PLAN.md`/`CHANGELOG.md` and, for breaking changes, a schema bump.
- Do not add dependencies without checking `AGENTS.md` §2.
- Database bundles must clear the license and size gates in `AGENTS.md` §6.

## Documentation

- User docs are bilingual (`docs/*.md` + `docs/*.zh.md`). Keep both in sync.
- Build with `pixi run -e docs mkdocs build --strict`; anchor/link warnings are
  treated as defects.
- Prefer generated facts (database tables, schema lists, counts) over hand-copied
  values; update the generator/tests when a fact changes.
- New public JSON documents should be registered in `gapit schema` where they are
  intended as agent contracts.

## Commits and releases

Use Conventional Commits (`feat:`, `fix:`, `test:`, `docs:`, `refactor:`,
`chore(release):`). Never commit unless the maintainer asks. Release tags run the
gates before publishing; keep `CHANGELOG.md` and `pyproject.toml` versions in
lockstep.

Report security issues through `SECURITY.md`, not public issues.
