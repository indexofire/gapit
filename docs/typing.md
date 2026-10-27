# Typing schemes

A typing document is the DATABASE-SIDE scoring spec that turns screening hits into phenotype
calls — serotypes, serogroups, pathotypes. Install one with `gapit db build --typing FILE`
and the database can designate any sample you screen against it. This page is the full
`gapit.typing/2` reference plus a cookbook: six researched designation schemes encoded as
validated, self-tested example documents.

- For the database-side mechanics (`--typing`, validation, manifests), see
  [Databases](./databases.md#typingjson-declarative-phenotype-scoring).
- For the designation output document (`gapit.typing_result/1`), see [Outputs](./outputs.md).
- The JSON Schema of a typing document prints with `gapit schema typing`.

## The two-stage designation workflow

Gene-kind designation is a two-stage pipeline (the rightsholder design): **`gapit screen`
detects genes, `gapit typing` designates from the results.**

```console
$ gapit screen contigs.fa --db ecoli_dec --output result.tsv --quiet
$ gapit typing result.tsv --quiet
FILE	SCHEME	PHENOTYPE	GENES	CONFIDENCE	SCORE	RUNNER_UP	NOTES
contigs.fa	gb4789_6	EHEC	escV;stx2a;uidA	high	1.0000	EAEC (1.0000)	...
contigs.fa	risk_monitoring	EHEC	escV;stx2a;uidA	high	1.0000	EAEC (1.0000)	...
```

Stage 1 is pure gene detection: a typed gene database screens byte-identically to an
untyped one — the frozen 15-column abricate table, in every format. Stage 2 reads one or
more gapit/abricate screen result tables (TSV or CSV), resolves the database from the
rows' shared `DATABASE` column, folds every `(FILE, GENE)` to its best row by
`(%IDENTITY, %COVERAGE)` (first row wins ties — the same per-gene fold the report hits
got), and evaluates the database's `typing.json` over each FILE's genes.

The two stages also compose straight through a pipe — the canonical one-liner:

```console
$ gapit screen 1.fna --db ecoli_dec --quiet | gapit typing --quiet
```

`gapit typing` reads its table from stdin whenever stdin is not a terminal (piped or
redirected); a bare invocation at a terminal keeps printing help. The explicit `-` marker
— `gapit typing -` — reads stdin even with a terminal attached, until EOF; mixing `-`
with file arguments is a usage error (v1 takes either stdin or files, never both). A
piped table's JSON `source` renders as `["-"]`.

Usage: `gapit typing RESULT.tsv [RESULT2.tsv ...] [-D datadir] [-f tsv|json|md] [-o FILE]
[-q]`, or piped as above (`-` for explicit stdin). The default TSV carries one row per
(FILE, scheme) with the columns `FILE`, `SCHEME`,
`PHENOTYPE`, `GENES`, `CONFIDENCE`, `SCORE`, `RUNNER_UP`, `NOTES` — ambiguous calls render
`-` for the phenotype (and runner-up) and carry the candidate pair in NOTES; `GENES`
follows the phenotype with the FILE's present gene names (sorted, `;`-joined), the same
list on each of the FILE's scheme rows — the hits that drove the call, in the open. `json` emits
`gapit.typing_result/1` (full score breakdowns plus the same per-file `genes` list;
`gapit schema typing_result`); `md` renders
the same eight columns under a frontmatter naming the db and source tables. Typed errors:
mixed `DATABASE` values are `DATABASE_MISMATCH`, a table with no data rows `TYPING_NO_DATA`
(the database is only knowable from rows), a database without a `typing.json`
`TYPING_NO_SCHEME`, and a cluster database `TYPING_CLUSTER_DB` — cluster typing never goes
through the command.

**The cluster exception.** A cluster database's designation stays integrated in
`gapit screen`: its unit of typing is the best *locus* call, which only the minimap2
cluster engine can compute — there is no re-readable screen table for it. Screening a
typed cluster database annotates the cluster document directly — a single-scheme db
annotates `gapit.cluster/1` (the typed-only `PHENOTYPE` TSV column, `best.phenotype` +
`best.phenotype_detail` in JSON/MD), and a multi-scheme db renders `gapit.cluster/2`
(one call per scheme in a `phenotypes` map, per-scheme winner-evidence TSV columns).

## One document, two engines

A `gapit.typing/2` document is a list of named **schemes**; each scheme is a list of
**rules** plus the decision thresholds applied to them:

```json
{
  "schema": "gapit.typing/2",
  "schemes": [
    {
      "name": "serogroup",
      "rules": [ /* one rule per candidate phenotype */ ],
      "cutoff": 0.9,
      "ambiguity_margin": 0.05,
      "fallback": "NG"
    }
  ]
}
```

Which engine evaluates the document depends on the database kind:

| Database kind | Input | Engine | Where the calls appear |
|---|---|---|---|
| **gene** | marker FASTA | blastn gene path | `gapit typing` over the screen results: `gapit.typing_result/1` (JSON) and the eight-column TSV/MD table |
| **cluster** (one scheme) | GBK/GFF loci | minimap2 locus path | integrated in `gapit screen`: `best.phenotype` + `best.phenotype_detail` of `gapit.cluster/1`, plus the typed-only `PHENOTYPE` TSV column |
| **cluster** (multi-scheme) | GBK/GFF loci | minimap2 locus path | integrated in `gapit screen`: `files[].phenotypes` of `gapit.cluster/2` (one call per scheme, compose rendered from its ingredients), plus per-scheme winner-evidence TSV columns |

Two hard constraints by kind:

- A **gene database** evaluates `weighted_genes` and `exact_set` rules only (no loci exist to
  feed `cluster_match`/`learned_linear`).
- A **cluster database** renders one of two contracts, picked by scheme count: a document
  with exactly **one** scheme (a `gapit.typing/1` document included) renders the
  byte-frozen `gapit.cluster/1` — one phenotype slot per file, `best.phenotype` — while a
  document with **more than one** scheme renders `gapit.cluster/2`: every scheme is called
  (schemes evaluate even without a best-locus call, so fallbacks always render), each file
  carries a `phenotypes` map keyed by scheme name in document order, `best` drops its
  phenotype fields, and the TSV/CSV columns derive from the document: each compose scheme's
  value column renders first (`SEROTYPE`), then every non-compose scheme carries its call
  column plus its winning locus's evidence (`<SCHEME>_LOCUS/_COVERAGE/_IDENTITY/_PRESENT/
  _PARTIAL/_MISSING`; the `o_group` evidence prefix shortens to `O`) — `-` cells where no
  winner exists. `/2` builds guarantee globally-unique
  gene ids: a `/gene` name shared across loci is resolved at build time to its
  `/locus_tag`, or to `locus:id@start` when no free locus_tag exists. A single
  compose-only scheme is unparseable (its placeholders must reference a non-composing
  sibling), so scheme count alone picks the contract.

## Rule and scheme primitives

Every rule is discriminated by its `model`. Full details live in `gapit schema typing`;
this is the working reference.

### Rule models

| Model | Scores | Key fields |
|---|---|---|
| `weighted_genes` | gene-presence sum, normalized to [0, 1] by the positive-weight total | `weights`, `negative` (count against the score when present), `identity_floor` (required), `coverage_floor` (optional), `require_any` (any-of gate that zeroes the rule when unsatisfied) |
| `exact_set` | 1.0 iff every `requires` gene is present AND at least one `requires_any` gene is present (when that set is non-empty) AND every `excludes` gene is absent, else 0.0 | `requires`, `requires_any` (any-of gate), `excludes`, floors defaulting to 90/90 |
| `cluster_match` | weighted coverage + identity + key-genes components over the best covered locus, each component floored | `coverage{weight,floor}`, `identity{weight,floor}`, `key_genes{weight,genes}` |
| `learned_linear` | sigmoid over weighted named features + bias | `features` (`gene:<id>:<cov\|ident\|present>`, `cluster:<locus>:<coverage\|identity>`), `weights`, `bias`, `trained_on` |

A gene "counts" as present for `weighted_genes`/`exact_set` when its best surviving hit (or
gene call, on the cluster path) clears the rule's floors. Keep `coverage + identity` weights
of a `cluster_match` rule below the cutoff: key genes should decide, so a fully covered but
key-gene-less locus cannot steal the call.

### Scheme fields

| Field | Meaning |
|---|---|
| `name` | `[a-z0-9_]+`, unique in the document; the designation output keys on it |
| `cutoff` | top rule must reach it to call at all, else the fallback |
| `ambiguity_margin` | top − second must reach it, else the call is `ambiguous` (null phenotype, top two listed) |
| `fallback` | low-confidence output string (`NT`, `NG`, `OUT`, `KUT`, …) |
| `control_gene` | gate: when this gene is absent, the scheme outputs its fallback with the note `control gene absent` |
| `unique_group` + `mixed_phenotype` | when more than one gene of a group is present, call `mixed_phenotype` with the pair in `ambiguous` (mixed-infection signal) |
| `compose` | `"{o_group}:{k_group}"`: a rule-less scheme rendered from sibling schemes' calls after they evaluate — fallback strings flow through, an ambiguous ingredient yields a null composition carrying that ingredient's variants |
| `notes` (rules) | free-text strings surfaced verbatim on the winning call — caveats, citations, flags |

One tie exception: a 1.0 tie that involves a satisfied `exact_set` resolves by **declaration
order**, not the margin (marker tables are mutually exclusive by construction, so a tie means
two rows matched and the earlier-declared one wins). Non-exact ties stay ambiguous.

## The cookbook

Six example documents live in
[`tests/data/typing/schemes/`](https://github.com/indexofire/gapit/tree/main/tests/data/typing/schemes),
one marker FASTA or GBK loci file each. Every document validates against the current schema,
and `tests/test_typing_schemes.py` / `tests/test_typing_schemes_cluster.py` build each fixture
database and assert the phenotype calls below. **All sequences are synthetic** — they encode
the scheme logic, not curated content; real curation is future work per provider licensing.
Two curated exceptions sit beside them: the bundled `lm_doumith` panel (real public-domain
INSDC sequences) and the bundled `cje_cps` capsule panel (amplicons extracted from public-domain
INSDC accessions), each with its own test module and described first below.

### C. jejuni penner capsule types — the bundled cje_cps panel (gene db)

`cje_cps` (37 records) is the capsule half of the updated Penner multiplex PCR, Poly et al.
2015 (PLoS ONE 10(12):e0144349, CC0): 36 amplicons — the exact PCR product of each published
primer pair, sliced from the accession the paper names — plus `lpxA` as the species
`control_gene`. 47 Penner serotypes collapse into 35 capsule types, so the scheme calls types
and complexes:

| Call | requires | notes |
|---|---|---|
| `HS4 complex (HS4/13/43/50/62/65)` | `HS4A` | the six serotypes share the locus; individual serotype needs antisera |
| `HS4/13/64 (CG8486)` | `HS4B` | CG8486 reference; `HS4B` also appears inside HS16/HS64 |
| `HS16` / `HS64` | `HS4A`+`HS4B`+`HS52` / `HS4A`+`HS4B` | the HS52 band is the only separator |
| `HS5/31` | `HS5/31`+`HS45` | shared CPS locus; HS45 band is an unexpected PCR product |
| `HS58` | `HS15`+`HS58` | HS15 band is shared with HS15 |
| `HS32` / `HS60` | `HS32`+`HS45` / `HS45`+`HS60` | HS45 band plus the locus-specific amplicon |
| `HS6/7`, `HS8/17`, `HS23/36`, `HS33/35` | one amplicon each | complex calls by design |
| every other type | its own amplicon | e.g. `HS2`, `HS3`, `HS10`, `HS41`, `HS45` |

Rules are `exact_set` (score 1.0 iff exactly the listed amplicons are present, with the rest
excluded), identity floor 95 / coverage floor 90 — the HS5/31 rule sits at 90 because its HS45
band is divergent. `lpxA` gates the scheme: without it the call is the fallback with the note
`control gene absent`. Fallback `untypeable (atypical CPS profile or non-jejuni)`.

Every caveat the paper documents lives in the winning rule's `notes` (the declared
unexpected amplicons, the HS1 607-vs-610 bp length difference, the complex groupings), so a
call carries its own qualification instead of requiring the paper to interpret it.

### Listeria serogroups — the bundled lm_doumith panel (gene db)

`lm_doumith` is a bundled database: the Doumith markers extracted from public-domain
INSDC records and the `doumith_serogroup` scheme (`gapit.typing/2`) — the same five-marker
multiplex-PCR semantics as LisSero, with `prs` as the genus-level `control_gene` and 95/95
floors on every rule (Doumith et al. 2004, J Clin Microbiol 42:3819):

| Call | requires | excludes |
|---|---|---|
| `IVb-v` | `ORF2819`, `ORF2110`, `lmo0737` | `lmo1118` |
| `IVb` (4b, 4d or 4e) | `ORF2819`, `ORF2110` | `lmo0737`, `lmo1118` |
| `IIb` (1/2b, 3b or 7) | `ORF2819` | `lmo0737`, `lmo1118`, `ORF2110` |
| `IIc` (1/2c or 3c) | `lmo0737`, `lmo1118` | `ORF2819`, `ORF2110` |
| `IIa` (1/2a or 3a) | `lmo0737` | `lmo1118`, `ORF2819`, `ORF2110` |

`IVb-v` — the 4b variant carrying the lineage-II `lmo0737` cassette (Huang 2011 "unusual
4b"; FDA 4bV; ST382/ST554 emergent clones) — is declared **before** `IVb` so a shared
ORF2819+ORF2110 profile carrying `lmo0737` resolves to the variant. Uppercase gene ids
(`ORF2819`, `ORF2110`) travel verbatim through build, screen, and typing. Two-stage use,
zero network (the bundle materializes on first screen):

```console
$ gapit screen isolates.fa --db lm_doumith --output lm.tsv --quiet
$ gapit typing lm.tsv --quiet
FILE	SCHEME	PHENOTYPE	GENES	CONFIDENCE	SCORE	RUNNER_UP	NOTES
lm_doumith_iia.fa	doumith_serogroup	IIa	lmo0737;prs	high	1.0000	IVb-v (0.0000)	...
lm_doumith_ivb.fa	doumith_serogroup	IVb	ORF2110;ORF2819;prs	high	1.0000	IVb-v (0.0000)	...
```

Fallback `untypeable (4a/4c, atypical profile, or non-Lm Listeria)`; without `prs` the
control gate fires (`control gene absent`). Known limits ride the calls as notes: 4b/4d/4e
share every marker, EGD-e types as IIc (documented in Doumith 2004), and HGT of the ORF
markers can fake a IIb profile — see the [lm_doumith section](databases.md) in the
databases guide.

### Listeria serogroups — the Doumith table (gene db)

`doumith.json` + `doumith.fa`: the complete five-marker multiplex-PCR table from Doumith et
al. 2004 (J Clin Microbiol 42:3819) — `prs` as the scheme `control_gene`, and one `exact_set`
rule per serogroup over `lmo0737`, `lmo1118`, `orf2819`, `orf2110`:

| Serogroup | requires | excludes |
|---|---|---|
| `1/2a-3a` | `lmo0737` | `lmo1118`, `orf2819`, `orf2110` |
| `1/2c-3c` | `lmo0737`, `lmo1118` | `orf2819`, `orf2110` |
| `1/2b-3b-7` | `orf2819` | `lmo0737`, `lmo1118`, `orf2110` |
| `4b-4d-4e` | `orf2819`, `orf2110` | `lmo0737`, `lmo1118` |

Fallback `NT` (4a/4c and atypical profiles). The 4b rule carries the caveat notes, including
the 4b* (IVb-v1) horizontal-gene-transfer trap: 4b* strains acquired the lineage-II
`lmo0737` cassette, match no plain pattern, and fall back to `NT` — the test suite pins
exactly that call (J Clin Microbiol 2022, PMID 36472431).

### Shigella / EIEC — ShigaTyper semantics (gene db)

`shigella.json` + `shigella.fa`: a skeleton of ShigaTyper (Wu et al. 2019, Appl Environ
Microbiol 85:e00165-19) with `ipaH_c` as the control gene and the wzx variants wired as one
`unique_group` (`mixed_phenotype: "mixed Shigella serotypes"` — ShigaTyper's "multiple wzx"
checkpoint becomes a mixed call). Highlights:

- **S. sonnei form I** (`exact_set` over `Ss_wzx`/`Ss_wzy`) is declared **before** the
  **form II** rule (`weighted_genes` on `Ss_methylase`): a form I sample whose antigen is
  phase-variable satisfies both at 1.0, and the exact-tie exception resolves to form I.
- **S. dysenteriae 1** requires `Sd1_wzx` + `Sd1_rfp`, with a note on the rfp-negative
  variant needing its own rule in a full panel.
- **S. flexneri** serotype conversions as exact sets over `gtr`/`Oac` (`2a`, `2b`, `3a`)
  plus a `Y/novel` catch-all (`requires` the Sf base, `excludes` the listed conversion
  genes) — with a full panel that rule splits into Y and ShigaTyper's "novel serotype".
- **EIEC** approximates ShigaTyper's checkpoint 3: `exact_set` requires `EclacY` and
  `excludes` the exemption genes `Sb9_wzx`/`Sb15_wzx`; `cadA` and ipaB plasmid-coverage
  nuances are documented as unmodelled in the note.

### N. meningitidis serogroups — meningotype panel (gene db)

`meningotype.json` + `meningotype.fa`: `ctrA` control gene, one `weighted_genes` rule per
serogroup over the Mothershed et al. 2004 (J Clin Microbiol 42:320–328) real-time panel
(`sacB`/A, `synD`/B, `synE`/C, `xcbB`/X, `synF`/Y, as named by the meningotype tool) — and
W/Y through the **EX7E allele probes** `synG_EX7E_P` / `synG_EX7E_G` at `identity_floor`
99.5 (see below). Fallback `NG`.

### V. parahaemolyticus O and K — the Kaptive pattern (cluster db)

`vp_ok.json` + `vp_ok.gbk`: OAgc/CPSgc-like loci after the Kaptive databases of
van der Graaf-van Bloois et al. 2023 (Microb Genom 9:mgen001007; 16 O and 71 K serotypes).
Three schemes:

- `o_group` — `cluster_match` rules (coverage 0.55 / identity 0.25 / key genes 0.2, floors
  95/95) with fallback `OUT` (O untypeable). The `O3/O13` rule encodes the combined label:
  the O3 and O13 loci contain the same genes and cannot be distinguished by gene content.
- `k_group` — the same shape over the K loci, fallback `KUT` (K untypeable).
- `serotype` — the `"{o_group}:{k_group}"` compose scheme documenting the panel output.

The vp_ok cookbook document installs as **one merged database** rendering
`gapit.cluster/2`: every scheme is called per file (`o_group`, `k_group`, then the
`serotype` compose rendered from its ingredients — `O3/O13:K6` in one scan), and its
TSV row puts the composed serotype FIRST and then gives each scheme its call plus the
winning locus's own evidence (coverage, identity, gene counts, missing ids):
`FILE SEROTYPE O_GROUP O_LOCUS O_COVERAGE ... K_GROUP K_LOCUS K_COVERAGE ...`. The
historical `/1` install —
**two databases** (`vp_o`, `vp_k`) built from the same GBK with one scheme each,
exactly how Kaptive ships separate O and K databases — remains fully valid: a
single-scheme document keeps rendering byte-frozen `gapit.cluster/1`, and the split
is still what you want when the O and K calls belong to different pipelines. The
test suite asserts the merged /2 goldens (`cluster2_vp_o3k6`, `cluster2_vp_none`)
and the split /1 calls (`O3/O13`, `K6`) side by side.

### V. cholerae O1/O139 and Ogawa/Inaba (cluster db)

`cholerae.json` + `cholerae.gbk`: loci `wbe` (O1: `wbeV`, `wbeW`, `wbeT`) and `wbm` (O139:
`wbmV`, `wbmW`, `wbfZ`) plus a non-typeable `wbc` locus, with two schemes:

- `serogroup` — `cluster_match` O1 / O139 with fallback `non-O1/non-O139`. The O139 rule's
  note records the **wbfZ trap**: the `wb*` cassette is exchanged through the conserved
  `gmhD`/`rjg` junctions, so environmental serogroups can carry O139-like cassettes and trip
  single-gene `wbfZ` assays (Sozhamannan et al. 1999, Infect Immun 67:6215) — the rule keys
  on the whole `wbm` locus instead.
- `subserotype` — Ogawa/Inaba through the wbeT allele probe: **ogawa** = `weighted_genes`
  on `wbeT` at `identity_floor` 99.9 (the determinant is single mutations in the wbeT
  methyltransferase; Stroeher et al. 1992, PNAS 89:2566; BMC Microbiol 2013, 13:173);
  **inaba** = `wbeV` positive with `wbeT: -2.0` in `negative` (the Inaba definition:
  wbeT not clearing the allele floor), with the Hikojima note — rare, unstable, no stable
  reference, not determinable. A mutated-wbeT sample calls `O1` on the serogroup database
  and `inaba` on the subserotype database at the same time (the `cluster_vc_inaba` golden).

### Diarrheagenic E. coli — the dual-scheme designation (gene db)

`dec.json` + `dec.fa`: the DEC designation as **two schemes over one gene database** — the
typing/2 multi-scheme showcase. Both schemes use exact semantics (cutoff 1.0, margin 0.0),
`uidA` as the `control_gene`, fallback `non-DEC`, and one severity-ordered rule ladder
(EHEC > STEC/EPEC > ETEC > EIEC > EAEC — hybrids surface as the runner_up, decided by
declaration order):

| Pathotype | requires | requires_any | excludes |
|---|---|---|---|
| `EHEC` | `escV` | `stx1a`/`stx1b`/`stx2a`/`stx2b` | |
| `STEC` | | `stx1a`/`stx1b`/`stx2a`/`stx2b` | `escV` |
| `EPEC` | `escV`, `bfpB` | | |
| `EPEC_atypical` | `escV` | | `bfpB`, every stx subunit |
| `ETEC` | | `lt`/`sth`/`stp` | |
| `EIEC` | `invE` | | |
| `EAEC` (`gb4789_6`) | | `aggR`/`pic`/`astA` | |
| `EAEC` (`risk_monitoring`) | `aggR` | | |

- **`gb4789_6`** encodes the GB 4789.6-2016 panel semantics: EAEC is *any of* aggR/pic/astA.
- **`risk_monitoring`** encodes the current food-safety risk-monitoring variant
  (最新食品安全风险监测方案): aggR is mandatory — pic/astA alone are insufficient. A
  pic+astA isolate without aggR is `EAEC` under gb4789_6 and `non-DEC` under
  risk_monitoring — the headline divergence between the two schemes, locked by the
  `ecoli_dec_typing` golden (both calls in one `gapit.typing_result/1` document; the
  bundled database carries the same document).
- `EPEC_atypical` is declared *after* EHEC/STEC so stx-positive isolates never land there;
  escV+/bfpB−/stx− isolates fall through the typical-EPEC rule onto it.
- The fixture mirrors the rightsholder database's duplicated records (pic ×2, sth ×3 there;
  one duplicate each in the fixture): duplicate records collapse by gene name — blastn's
  `-culling_limit 1` keeps the best subject per query span and the engine folds one gene
  call per name.

## The allele probe trick (and its successor)

Two cookbook schemes discriminate alleles — not genes — using presence/absence primitives:

- **EX7E (meningotype)**: the synG EX7E motif peptide reads P in serogroup W, G in Y, S in
  dual W/Y (meningotype's `menwy` check). The fixture carries the two allele sequences as
  separate probe records at `identity_floor` 99.5: only the exact allele clears the floor,
  because a single SNP in a ~120 bp probe already drops below 99.5% identity.
- **wbeT (cholerae)**: Ogawa vs Inaba is single mutations in wbeT, so ogawa requires wbeT
  at 99.9 and inaba is the negative rule — an inactivating mutation (98% identity) still
  maps, still counts as a present gene for the locus rules, but fails the allele floor.

Two practical constraints:

- The probe records must **stand alone** in the database. gapit runs blastn with abricate's
  `-culling_limit 1` (parity is a feature), and a long record containing the probe region
  would shadow the probe's alignment — so the meningotype fixture has no full-length synG
  record next to the probes.
- An allele *between* the probes (the dual W/Y S motif, or a novel wbeT missense) clears
  neither floor and lands on the fallback — allele-level ambiguity, surfaced as a low
  confidence call.

When a future `allele_match` primitive lands (pattern + allele table, judged on the aligned
query sequence the way meningotype translates its EX7E window), both tricks supersede to a
single rule with explicit alleles.

## Calibration and the learned_linear outlook

Tuning thresholds against labeled assemblies is what
[`scripts/cluster_calibration.py`](https://github.com/indexofire/gapit/blob/main/scripts/cluster_calibration.py)
is for (developer tool): point it at a typed cluster database and a label TSV and it prints
the per-expected-phenotype score distribution, the called×expected agreement matrix, and the
divergence list:

```console
$ python scripts/cluster_calibration.py mycps --datadir ./db labels.tsv
```

That report is the loop a future `learned_linear` training flow consumes: the feature
vectors (`gene:<id>:cov|ident|present`, `cluster:<locus>:coverage|identity`) are exactly the
calibration columns, so a trained model drops into the same scheme as one more rule.

## Building and typing a cookbook scheme

```console
$ gapit db build doumith doumith.fa --typing doumith.json
{"db":"doumith","records":5,"dbtype":"nucl","destination":".../doumith"}

$ gapit screen lm_4b.fa --db doumith --output lm_4b.tsv --quiet
$ gapit typing lm_4b.tsv --quiet
FILE	SCHEME	PHENOTYPE	GENES	CONFIDENCE	SCORE	RUNNER_UP	NOTES
lm_4b.fa	doumith	4b-4d-4e	orf2110;orf2819;prs	high	1.0000	1/2a-3a (0.0000)	...

$ gapit screen vc_inaba.fa --db vc_subserotype --quiet
FILE      BEST_LOCUS  TYPE  PHENOTYPE  COVERAGE  IDENTITY  PRESENT  PARTIAL  MISSING_IDS
vc_inaba.fa  wbe      O1    inaba      100.00    99.44     3        0        -
```

The second block is the cluster exception: a typed cluster database designates inside
`gapit screen` (no second stage — there is no screen table to re-read).

The cluster O/K pattern: the merged /2 install (one scan, every scheme) and the
historical /1 split (one database per scheme):

```console
$ gapit db build vp_ok vp_ok.gbk --typing vp_ok.json   # ONE merged db → gapit.cluster/2
$ gapit screen sample.fa --db vp_ok --quiet
FILE      SEROTYPE   O_GROUP  O_LOCUS     O_COVERAGE  ...  K_GROUP  K_LOCUS  K_COVERAGE  ...
sample.fa O3/O13:K6  O3/O13   OL3_or_O13  100.00      ...  K6       KL6     100.00      ...

$ gapit db build vp_o vp_ok.gbk --typing vp_o_group.json   # /1 split (one scheme each)
$ gapit db build vp_k vp_ok.gbk --typing vp_k_group.json
$ gapit screen sample.fa --db vp_o --quiet   # → O3/O13
$ gapit screen sample.fa --db vp_k --quiet   # → K6
```
