# PLAN.md — gaita development roadmap

> Phase-gated. Each phase ends green on `pixi run lint && pixi run typecheck && pixi run test`.
> Behavior targets come from `SPEC.md`; engineering rules from `AGENTS.md`.

## Phase 0 — Scaffolding

**Goal**: runnable, gated, empty skeleton.

- pixi deps: `python =3.14`, `blast >=2.7`, `any2fasta`; PyPI: `typer`, `pydantic >=2`, `rich`;
  dev: `ruff`, `basedpyright`, `pytest`.
- Package metadata: `requires-python = ">=3.11"` (develop on 3.14, support 3.11+; CI matrix
  3.11 / 3.13 / 3.14).
- pixi tasks: `lint`, `fmt`, `typecheck`, `test`, `gaita`. src-layout package with `py.typed`.
- `gaita --version` prints and exits 0.

**Done when**: all gates pass on the skeleton; `gaita --version` works.

## Phase 1 — Database layer (`fasta.py`, `db.py`, `config.py`)

**Goal**: gaita understands an abricate datadir.

- Streaming FASTA reader (plain + gz/bz2 via `any2fasta` fallback), strict header model.
- DB discovery (`<datadir>/<name>/sequences`), `~~~` header parsing with upstream fallback rules
  (SPEC §4 step 5), `mol_type` heuristic (SPEC §2), `makeblastdb` wrapper (argv list).
- CLI: `gaita list` (DATABASE/SEQUENCES/DBTYPE/DATE table + `--json`), `gaita setupdb`.

**Done when**: unit tests for header parsing/mol_type; `gaita list` against a fixture datadir
matches `abricate --list` output.

## Phase 2 — Screening core (`blast.py`, `hits.py`, `minimap.py`, `report.py`)

**Goal**: the parity engine. TDD: tests written from SPEC.md first.

- any2fasta→blastn pipeline (exact flags, SPEC §3), 15-field row parser (typed at the boundary).
- Hit processing in SPEC order: strand swap → `(qseqid,qstart,qend)` dedup → unrounded coverage
  filter → field fallbacks → product cleanup. **No interval merging.**
- `minimap` replicated with exact int arithmetic (SPEC §4), incl. the 1-based/box-0 quirk.
- `Report` pydantic model as the canonical in-memory result; stable sort SEQUENCE/START.

**Done when**: unit tests cover every formula and quirk in SPEC §4/§7 (incl. the 79.996% filter
case, partial `~~~` headers, minus-strand dedup-collapse).

## Phase 3 — TSV/CSV output + parity harness

**Goal**: byte-compatible drop-in replacement.

- `formats/tsv.py`: header rules, `--noheader`, `--nopath`, `--csv`, per-file buffering with
  preserved row order.
- `pixi run parity`: run real abricate (conda) and gaita over a small committed genome corpus
  against ≥2 DBs (ncbi + card); diff gene calls byte-for-byte (TSV) modulo the FILE column.
- CLI flags complete for screening: `--db --datadir --minid --mincov --threads --fofn --quiet
  --csv --noheader --nopath --debug`.

**Done when**: parity diff is empty on the corpus. **This is the v1.0 release gate.**

## Phase 4 — Agent outputs (`formats/json.py`, `formats/md.py`, `errors.py`)

**Goal**: the reason gaita exists.

- `--format json|md` on the screening path; top-level `"schema": "gaita.report/1"`,
  snake_case keys, explicit units (`identity_pct`, `coverage_pct`).
- Markdown: YAML frontmatter (version, db, params, ISO-8601 UTC) + stable tables.
- Typed errors → JSON envelope `gaita.error/1` on stderr + documented exit codes (SPEC §1).
- Self-description: `gaita schema report|summary|error` (prints JSON Schema from the pydantic
  models), `gaita list --json`, `gaita --version --json`.

**Done when**: golden files for JSON/MD; `gaita schema` output validates against the models;
error paths produce the envelope (tests force each exit code).

## Phase 5 — Summary mode (`summary.py`)

**Goal**: `gaita summary` matrix, parity with `abricate --summary`.

- Dutch mode (single multi-FILE report), union-of-genes columns, `;`-joined cells, `.` absent,
  `NUM_FOUND` distinct-gene count, `--identity`, zero-hit files included.
- Matrix also emitted as JSON (`gaita.summary/1`) and Markdown `[gaita-extension]`.

**Done when**: matrix parity against `abricate --summary` on the Phase-3 corpus; golden JSON/MD.

## Phase 6 — DB acquisition (`gaita db fetch`) — post-1.0

- Reimplement `abricate-get_db` per DB (ncbi, card, resfinder, argannot, plasmidfinder, megares,
  ecoh, vfdb, ecoli_vf, bacmet2, victors, upec_expec_vf) with the documented transforms
  (SPEC §8). Until then, gaita consumes abricate-built datadirs.
- Optional: ship a pixi-packaged snapshot of the bundled DBs.

## Phase 7 — Hardening & distribution — post-1.0

- Conda/pixi package recipe; shell completions; `--debug` parity of stderr diagnostics.
- Performance pass (multi-file parallelism across inputs; threads >1 determinism check vs
  abricate behavior).
- Optional: MCP server exposing `screen`/`summary`/`schema` for agent runtimes.

## Milestones

| Milestone | Contents | Gate |
|---|---|---|
| **M1 (alpha)** | Phases 0–2 | unit suite green; screens a fixture genome |
| **M2 (beta)** | Phase 3 | parity diff empty |
| **v1.0** | Phases 4–5 | agent contract frozen as `*/1`; full gates green |

## Open decisions (resolve before the phase that depends on them)

1. **Protein DBs (blastx path)**: reproduce upstream's ignored-`--minid` quirk (current plan:
   yes, with stderr note) — decide before Phase 2 closes.
2. **JSON field naming for COVERAGE_MAP / GAPS**: keep abricate strings (`coverage_map`,
   `gaps = "openings/gaps"`) vs structured objects — decide in Phase 4; v1 keeps the strings for
   1:1 mapping, structured views can land in `gaita.report/2`.
3. **CSV + summary format mixing**: gaita auto-detects separator per report file
   `[gaita-extension]`; confirm no parity test relies on the broken upstream behavior.
