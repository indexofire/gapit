# Databases

A gapit database is a directory of reference genes you screen contigs and reads against.
This page covers where databases live, how to get them, what gapit writes on disk, and how
to bring your own. For the full datadir layout, sharing, conda/Docker guidance, and the
database anatomy, see [Data directory](./datadir.md). For how screening uses them, see
[screen.md](./screen.md).

## Database naming

New database names describe the organism with a short code followed by a content suffix:

| Scope | Rule | Examples |
|---|---|---|
| Species | genus initial + first two species letters | `cje_cps` (*C. jejuni* capsule), `eco_vf` (*E. coli* virulence), `vch_serogroup` (*V. cholerae* serogroup), `vpa_ok` (*V. parahaemolyticus* O:K) |
| Genus | first three genus letters | `shi_serotype` (*Shigella*), `cam_cps` (*Campylobacter* genus-level capsule) |

Use lowercase ASCII, underscores, and a content suffix (`_cps`, `_vf`, `_serogroup`,
`_ok`, `_dec`). Existing published names (`cam`, `ecoh`, `ecoli_dec`, `lm_doumith`,
`vc_serogroup`, `vp_ok`, ...) are legacy identifiers: they stay stable for compatibility
and are not renamed by this convention.

## Where databases live

Every database sits in a **datadir**, one subdirectory per database. gapit resolves the
datadir in this order:

1. `--datadir` on the individual call
2. `$GAPIT_DATADIR`
3. `$XDG_DATA_HOME/gapit/db` (absolute `XDG_DATA_HOME` only)
4. `~/.local/share/gapit/db`

The location is resolved even before it exists. `gapit db list` and `gapit setupdb`
bootstrap the missing root automatically (so a fresh `pip install gapit` can be inspected
and set up immediately); `db fetch`/`db build` create it as part of their write path.
Commands that consume installed databases require it — except that `screen`/`typing`
against a **bundled** name materialize that bundle and create the root on first use. A
missing root on a strict read path is a typed error and exit 4:

```console
$ gapit db outdated --datadir /no/such/dir
{"schema":"gapit.error/1","code":"DATADIR_NOT_FOUND","message":"datadir does not exist: /no/such/dir","context":{"datadir":"/no/such/dir"}}
```

To keep databases inside an activated conda environment instead of your home directory,
opt in explicitly (recommended only when you need env-local data — the home default
survives environment rebuilds and is shared by every env):

```bash
export GAPIT_DATADIR="${CONDA_PREFIX}/db/gapit"
```

## Listing what's installed

`gapit db list` shows every known database and whether it is installed:

```console
$ gapit db list
NAME	DESCRIPTION	STATUS	DBTYPE	PROVIDER
ab_k	A. baumannii K locus (Kaptive)	available	nucl	Kaptive (Kenyon lab)
ab_o	A. baumannii OC locus — official keyword ab_o (Kaptive)	available	nucl	Kaptive (Kenyon lab)
argannot	ARG-ANNOT acquired resistance genes	available	nucl	IHU Méditerranée-Infection
bacmet2	BacMet2 experimentally confirmed biocide/resistance genes (protein)	available	prot	University of Gothenburg
card	CARD protein homolog resistance models	available	nucl	McMaster University
ecoli_kps	E. coli group 2+3 capsular polysaccharide loci (Kaptive)	available	nucl	Kaptive (Gladstone lab)
ecoli_vf	E. coli virulence factors (phac-nml)	available	nucl	PHAC-NML
kosc_k	K. oxytoca species complex K locus (Kaptive)	available	nucl	Kaptive (klebgenomics)
kosc_o	K. oxytoca species complex O locus (Kaptive)	available	nucl	Kaptive (klebgenomics)
kpsc_k	K. pneumoniae species complex K locus (Kaptive)	available	nucl	Kaptive (klebgenomics)
kpsc_o	K. pneumoniae species complex O locus (Kaptive)	available	nucl	Kaptive (klebgenomics)
megares	MEGARes antimicrobial resistance genes	available	nucl	MEG Lab
plasmidfinder	CGE PlasmidFinder replicons	available	nucl	DTU CGE
upec_expec_vf	UPEC/ExPEC virulence genes (FordeGenomics)	available	nucl	FordeGenomics
vfdb	VFDB virulence factors (set A, nucleotide)	available	nucl	USTC (VFDB)
victors	Victors virulence factors	available	nucl	University of Chicago
cje_cps	Campylobacter jejuni Penner capsule (CPS) multiplex PCR panel	bundled	nucl	gapit-curated (Poly 2015; CC0 + NCBI public domain)
ecoh	E. coli O and H antigens (srst2 EcOH)	bundled	nucl	Holt lab (srst2)
ecoli_dec	Diarrheagenic E. coli marker panel (GB 4789.6 + risk-monitoring designation)	bundled	nucl	gapit-curated (public-domain sources)
lm_doumith	Listeria monocytogenes serogrouping (Doumith 2004)	bundled	nucl	gapit-curated (public-domain INSDC sources)
ncbi	NCBI AMRFinderPlus (reference finder) curated AMR	bundled	nucl	NCBI
resfinder	CGE ResFinder acquired resistance genes	bundled	nucl	DTU CGE
vc_serogroup	Vibrio cholerae O1/O139 serogroup + ctxA toxigenicity markers; handed over from gside	bundled	nucl	gapit-curated (NCBI X59554.1 slice / KY660230.1 / X00171.1)
vp_ok	V. parahaemolyticus O and K antigen loci (OAgc/CPSgc, vpasero)	bundled	nucl	gapit-curated (vpasero)
```

NAME is the database name you pass to `--db`; PROVIDER names the upstream maintainer
organisation. STATUS reads `installed (N)` when `<datadir>/<name>/gapit-manifest.json`
exists, with N the record count, otherwise `available`. The last eight rows are the
wheel-shipped [bundled](#bundled-databases-install-time-ready) databases, alphabetical:
`bundled` before materialization, `installed (N)` after. Three of them (`ecoh`, `ncbi`,
`resfinder`) are also registry providers — such names render
exactly once, in the bundled section (never as a duplicate registry `available` row), and
`gapit db fetch <name>` remains their fresh-upstream update path; `vp_ok` is bundled-only
(no provider, no fetch). Databases installed
into the datadir outside the catalog also appear, after the bundled section and sorted by
name: `db build` products (gene and cluster kinds alike) and manifest-less directories
(abricate-style or `db install` bytes) render with PROVIDER `local`, the record count
from the manifest — or counted from the FASTA when no manifest exists — and DBTYPE from
the manifest, the BLAST index suffix, or abricate's letter heuristic, in that order.
On an interactive terminal the same rows render as a styled rich table; piped or
redirected output always stays the plain TSV above. For agents, `--json` emits a
`gapit.dblist/1` document (first three entries shown, output trimmed):

```console
$ gapit db list --json
{
  "schema": "gapit.dblist/1",
  "providers": [
    {
      "name": "argannot",
      "vendor": "IHU Méditerranée-Infection",
      "description": "ARG-ANNOT acquired resistance genes",
      "dbtype": "nucl",
      "installed": true,
      "kind": "gene",
      "records": 2224
    },
    {
      "name": "bacmet2",
      "vendor": "University of Gothenburg",
      "description": "BacMet2 experimentally confirmed biocide/resistance genes (protein)",
      "dbtype": "prot",
      "installed": true,
      "kind": "gene",
      "records": 746
    },
    {
      "name": "card",
      "vendor": "McMaster University",
      "description": "CARD protein homolog resistance models",
      "dbtype": "nucl",
      "installed": true,
      "kind": "gene",
      "records": 6059
    }
  ]
}
```

| Field | Type | Meaning |
|---|---|---|
| `schema` | string | Always `gapit.dblist/1` |
| `providers` | array | One entry per known database |
| `providers[].name` | string | Database name, as passed to `--db` / `gapit db fetch` |
| `providers[].vendor` | string | Upstream maintainer organisation (NCBI, DTU CGE, Kaptive (klebgenomics), ...) |
| `providers[].description` | string | Short content summary |
| `providers[].dbtype` | string | `nucl` (screened with blastn) or `prot` (screened with blastx) |
| `providers[].kind` | string | `gene` (per-gene BLAST screening) or `cluster` (minimap2 locus screening, [screen.md](./screen.md#cluster-databases-kind-cluster)) |
| `providers[].installed` | boolean | True when a manifest and readable `sequences` exist in the datadir |
| `providers[].indexed` | boolean | True when a BLAST index (`.nin` or `.pin`) exists; `setupdb` rebuilds a missing/stale one |
| `providers[].records` | integer | Record count, omitted when the database isn't installed |
| `providers[].license` | string | Upstream content license, omitted unless the provider or bundled snapshot pins one (card, vfdb, ecoli_vf, kaptive, bundled content) |
| `providers[].source` | string | `local` for a datadir-discovered database outside the catalog (`db build` / `db install`), `bundled` for a wheel-shipped database (see below); absent for registry entries |

`gapit db list` is the single listing surface: the database catalog above, with `--json`
returning the `gapit.dblist/1` document. (The former standalone listing command and its
schema were removed; abricate `--list` byte-parity is intentionally dropped for this
surface.) See [outputs.md](./outputs.md).

## Inspecting one database

`gapit db list` answers "which databases exist, and are they installed";
`gapit db info NAME` answers "what is this database, where does its content come from,
how is it meant to be used, and where does it mislead". It is a pure read over four
existing sources — the provider catalog, the wheel's `bundled.json`, the installed
manifest, and an optional `info.json` sidecar in the database directory — so nothing
has to be re-declared:

```console
$ gapit db info cje_cps
NAME	cje_cps
DESCRIPTION	Campylobacter jejuni Penner capsule (CPS) multiplex PCR panel
VENDOR	gapit-curated (Poly 2015; CC0 + NCBI public domain)
KIND	gene
DBTYPE	nucl
STATUS	installed
RECORDS	37
DESTINATION	/home/you/.local/share/gapit/db/cje_cps
SNAPSHOT	2026-10-27
TYPING_SCHEMA	gapit.typing/2
TOOL	gapit 0.6.4
FETCHED_AT	2026-10-27T04:19:06Z
MAKEBLASTDB_VERSION	makeblastdb: 2.17.0+
MINIMAP2_VERSION	2.31-r1302
SOURCE_URL	local
REFERENCE	Poly F, Serichantalergs O, Kuroiwa J, ... PLoS ONE 2015;10(12):e0144349.
REFERENCE_DOI	10.1371/journal.pone.0144349
REFERENCE_URL	https://doi.org/10.1371/journal.pone.0144349
REFERENCE_NOTE	Source of the 36 primer pairs and their expected amplicon sizes (Table 3) ...
REFERENCE	Klena JD, Parker CT, Knibb K, ... J Clin Microbiol 2004;42(12):5549-57.
REFERENCE_NOTE	Source of the lpxA positive-control primer set used as the panel's species gate.
CAVEAT	47 Penner serotypes collapse into 35 capsule types: ...
CAVEAT	Short amplicons are E-value bound, not --minid bound: ...
USAGE	gapit screen isolate.fa --db cje_cps --output cps.tsv   # stage 1: detect CPS amplicons
USAGE	gapit typing cps.tsv                                   # stage 2: penner_cps call
```

The table is a stable two-column `FIELD<TAB>VALUE` stream: list-valued fields
(`REFERENCE`, `CAVEAT`, `USAGE`, `SOURCE_URL`) repeat their key, and a reference's
`REFERENCE_DOI`/`REFERENCE_URL`/`REFERENCE_NOTE` rows follow its own `REFERENCE` row.
Fields a source does not carry are omitted, so the shape stays honest: a
bundled-but-unmaterialized database has no `TOOL`/`RECORDS` rows, and an uninstalled
provider has neither build provenance nor a record count.

`--json` emits the `gapit.dbinfo/1` document (`gapit schema dbinfo`) for agents:

```console
$ gapit db info vc_serogroup --json | head -3
{"schema":"gapit.dbinfo/1","name":"vc_serogroup","description":"Vibrio cholerae O1/O139 ...
```

`STATUS` mirrors `db list` — `installed`, `bundled` (ships in the wheel, materializes on
first use), or `available` (fetch-time download). An unknown name is exit 4
`DATABASE_NOT_FOUND`; a broken sidecar is exit 4 `DBINFO_MALFORMED`.

**Bringing your own literature.** A `db build` product can carry its own
`<datadir>/<NAME>/info.json` with the same three optional lists:

```json
{
  "references": [{"citation": "Smith 2026, in-house panel (unpublished)"}],
  "caveats": ["unpublished markers: confirm before reporting"],
  "usage": ["gapit screen x.fa --db labgenes"]
}
```

The sidecar wins over a bundled panel's declared guidance, so a datadir copy can
tighten or replace the wheel's notes without rebuilding anything.

## Bundled databases (install-time ready)

Most databases download from upstream at `db fetch` time because their licenses forbid
redistribution (see [Providers](#providers)). A small set of **bundled** databases ships
inside the gapit wheel instead — eight databases whose content provenance was audited
record-by-record (public-domain sources, or Apache-2.0 / BSD-3-Clause / MIT panels;
GPL and non-commercial content never rides the wheel). Nothing to download, nothing to
accept: a bundled database works out of the box.

| Name | Content | Snapshot | Typed |
|---|---|---|---|
| `cje_cps` | C. jejuni Penner capsule multiplex PCR panel (37 records: 36 CPS amplicons + `lpxA` control) | 2026-10-27 | `gapit.typing/2` (scheme `penner_cps`) |
| `ecoh` | E. coli O and H antigens (597 records, srst2 EcOH) | 2026-10-02 | — |
| `ecoli_dec` | Diarrheagenic E. coli marker panel (16 records, GB 4789.6 + risk-monitoring designation) | 2026-10-02 | `gapit.typing/2` (schemes `gb4789_6`, `risk_monitoring`) |
| `lm_doumith` | Listeria monocytogenes serogrouping (5 records, Doumith 2004 markers) | 2026-10-03 | `gapit.typing/2` (scheme `doumith_serogroup`) |
| `ncbi` | NCBI AMRFinderPlus curated AMR (8373 records) | 2026-10-02 | — |
| `resfinder` | CGE ResFinder acquired resistance genes (3206 records) | 2026-10-02 | — |
| `vc_serogroup` | V. cholerae O1/O139 serogroup + ctxA toxigenicity markers | 2026-10-19 | `gapit.typing/2` (scheme `serogroup`) |
| `vp_ok` | V. parahaemolyticus O and K antigen loci (195 loci, vpautils OAgc+CPSgc; cluster kind) | 2026-10-12 | `gapit.typing/2` (schemes `o_group`, `k_group`, `serotype`) |

**License provenance.** Every `ecoli_dec` record was re-sourced from primary
public-domain submissions (NCBI RefSeq/GenBank/DDBJ accessions, same alleles as the
rightsholder panel); two VFDB-derived records found in the original panel were removed
and replaced during the 2026-10 audit, and the byte-identical `estA4` (J03311.1) copy of
`estA3` (M18346.1) — one `sth` gene, one sequence — was dropped by the build pipeline's
exact-duplicate dedup in the 2026-10 rebuild, leaving 16 records. The `lm_doumith` panel was extracted the same
way from public-domain INSDC records (accessions in each header; reverse-complement
markers strand-corrected, `lmo0737` verified 100% against EGD-e). The three provider
snapshots (`ecoh`, `ncbi`, `resfinder`) passed the 2026-10
content-level audit across **all** records: NCBI AMRFinderPlus content is public
domain (US government work), the CGE ResFinder database is Apache-2.0, and srst2's EcOH
is BSD-3-Clause — all three permit
redistribution inside the MIT-licensed wheel. The `vp_ok` cluster snapshot is
MIT-licensed vpautils content (lab-curated loci, Bakta-annotated, via
[indexofire/vpautils](https://github.com/indexofire/vpautils)); it ships as the
resolved `sequences` + `features.json` pair a completed gapit.typing/2 build
produces. The `cje_cps` panel is assembled from the GenBank accessions the source paper
lists (public-domain INSDC records; the paper itself is CC0), so it clears the same gate.
A regression test pins that no `VF*` / `VFDB` / `ARO:` tags appear in any
of the eight bundles' headers.

**Snapshots are point-in-time.** Each snapshot's `sequences` file is byte-for-byte what
the upstream content produces on the snapshot date (a regression test proves the
byte-identity for `resfinder` against a pinned copy of the upstream archive, and for
`vp_ok` against a `gapit db build` over the real vpautils GenBank files when local
copies exist — `features.json` included, gene-id resolution included).
The wheel copy therefore drifts as upstream curation moves on: for the provider-backed
names (`ecoh`, `ncbi`, `resfinder`) `gapit db fetch <name>`
(alone, or with `--force` to overwrite a materialized copy) re-downloads the **latest**
upstream content and rebuilds the datadir database — the freshness path. `vp_ok` has no
provider: a refresh means a new resolved snapshot riding a wheel release.
`bundled.json`
records the snapshot date for reference. (`db outdated` stays manifest-based in v1:
it reports the materialization/fetch time of what is installed, not bundled-vs-upstream
drift — comparing the wheel snapshots against latest upstream is future work.)

**How materialization works.** The first `gapit screen ... --db <bundled-name>` against
a datadir that does not hold the database builds it there automatically — one stderr note
(`gapit: materializing bundled database ecoli_dec (16 records) into <datadir>`,
silenced by `--quiet`), then the standard gene-build pipeline: `records.jsonl`, the
`sequences` FASTA with gapit/v1 headers, the BLAST index, the `typing.json` copy when the
bundle ships one, and a manifest stamped `source: "bundled"`. That pipeline — shared with
`gapit db build` — drops exact `(gene, sequence)` duplicate records at build time (first
kept; a stderr note and the receipt's `duplicates_dropped` count the drops). Cluster
bundles (`vp_ok`) run the cluster pipeline instead: materialization rebuilds the locus
FASTA and `features.json` from the committed artifacts through the same cluster build
entry a `db build` uses — gene ids stay snapshot-resolved, the typing spec re-validates
against them, no re-resolution. Deterministic, zero
network. A missing datadir is bootstrapped by `gapit db list`, `gapit setupdb`,
`db fetch`/`db build`, and the first screen against a bundled name; `gapit setupdb`
materializes every bundled database alongside indexing; re-running either is a no-op
once the manifest exists. `db list` shows the database as `bundled` before that and
`installed (16)` after — on all four listing surfaces (TSV, rich table, `--json`, MCP
`db_list`).

### C. jejuni Penner capsule typing (`cje_cps`)

`cje_cps` is the *Campylobacter jejuni* capsule (Penner) panel: the 36 capsule-locus
amplicons of the updated multiplex PCR in Poly et al. 2015 (PLoS ONE 10(12):e0144349,
Table 3) plus an `lpxA` species control. Each record is the **amplicon** the published
primer pair produces — the exact sequence the PCR would amplify — extracted from the
GenBank accession the paper names for that primer set (the `acc` field carries the
accession and slice). 47 Penner serotypes collapse into 35 capsule types, so the panel
targets types and complexes rather than individual serotypes.

Designation is two-stage over a gene database (`gapit.typing/2`, scheme `penner_cps`,
`lpxA` as `control_gene`):

```console
$ gapit screen isolate.fa --db cje_cps --output cj.tsv --quiet
$ gapit typing cj.tsv --quiet
FILE	SCHEME	PHENOTYPE	GENES	CONFIDENCE	SCORE	RUNNER_UP	NOTES
isolate.fa	penner_cps	HS4 complex (HS4/13/43/50/62/65)	HS4A;lpxA	high	1.0000	HS16 (0.0000)	...
```

What the panel can and cannot resolve is encoded in the rules, not hidden:

- **Complexes are single calls** — `HS4 complex` (4/13/43/50/62/65), `HS6/7`, `HS8/17`,
  `HS23/36`, `HS33/35`, `HS5/31` share their CPS locus at the nucleotide level, so the
  call names the complex. An individual serotype inside a complex needs antisera or
  cgMLST.
- **Shared amplicons are profiles, not noise** — HS5/31 also amplifies the HS45 band,
  HS58 shares the HS15 band, HS32/HS60 produce an HS45 band, and HS16/HS64 are separated
  only by the HS52 band. Those multi-band profiles are declared as complete `exact_set`
  rules (HS16 = `HS4A+HS4B+HS52`, HS64 = `HS4A+HS4B`, HS5/31 = `HS5/31+HS45`,
  HS58 = `HS15+HS58`, HS32 = `HS32+HS45`, HS60 = `HS45+HS60`), each carrying the paper's
  caveat in its notes.
- **`lpxA` gates everything** — without a positive `lpxA` call the scheme reports the
  fallback with the note `control gene absent`, because the multiplex only means anything
  on a *C. jejuni* template.
- **One documented length deviation** — the HS1 primer pair amplifies 607 bp in the
  reference locus (BX545859.1) where the paper reports 610 bp; the record keeps the actual
  reference slice and its bundled note states the difference.

Fallback `untypeable (atypical CPS profile or non-jejuni)`. Coverage floor 90, identity
floor 95 (the HS5/31 rule sits at 90 because the HS45 band in that locus is a divergent
amplicon).

`gapit db info cje_cps` carries the panel's four literature references (Poly 2015 for the
panel, Klena 2004 for the `lpxA` control, Penner 1980 for the scheme, Pike 2013 for the
serotype distribution) and the seven interpretation caveats — the complex collapse, the
shared-amplicon profiles, the `lpxA`-is-not-speciation warning, the 62 bp HS2 E-value
bound, the HS1 607-vs-610 bp difference, the HGT caveat, and the paper's own
untypeable-rate observation.

### Listeria serogrouping (`lm_doumith`)

The typed bundle behind Listeria monocytogenes serogroup prediction: five markers —
`prs`, `lmo0737`, `lmo1118`, `ORF2819`, `ORF2110` — screened with 95/95 identity and
coverage floors and folded by the `doumith_serogroup` scheme into Doumith's multiplex-PCR
groups (Doumith et al. 2004, J Clin Microbiol 42:3819), extended with the Huang 2011
4b variant: `prs` gates the scheme as a genus-level control gene, `IIa` (1/2a or 3a),
`IIc` (1/2c or 3c), `IIb` (1/2b, 3b or 7), `IVb` (4b, 4d or 4e), and `IVb-v` (the
lmo0737-carrying 4b variant) — the same in-silico semantics as tools like LisSero.
Designation runs through the two-stage pipeline: `gapit screen -o table.tsv --db
lm_doumith`, then `gapit typing table.tsv`.

Known limitations, carried on the calls as notes where relevant:

- **4b/4d/4e are unresolvable at gene level** — the three serovars share every published
  molecular marker; resolve by cgMLST or antisera.
- **IVb-v is an emergent concern** — lmo0737 in an otherwise-4b profile (Huang 2011
  "unusual 4b"; ST382/ST554 clones) is declared before IVb so the variant wins.
- **EGD-e types as IIc** — the 1/2a reference genome carries lmo1118 and lands in IIc;
  documented in Doumith 2004, not a tool artifact.
- **Horizontal gene transfer cuts both ways** — the lineage-II lmo0737 cassette moves
  between lineages, and rare HGT of the ORF markers can produce false IIb-looking
  profiles; treat single-marker IIb calls in epidemiologically unlikely contexts with
  care.
- **`prs` is genus-level, not Lm-specific** — other Listeria species carry it; a prs+
  isolate with no serogroup markers falls back rather than being called non-Listeria.

The fallback call — `untypeable (4a/4c, atypical profile, or non-Lm Listeria)` — is the
honest answer for everything the five markers cannot resolve.

### vp_ok (V. parahaemolyticus O:K, MIT, typed on install)

`vp_ok` wraps the two *Vibrio parahaemolyticus* serotype references of
[indexofire/vpautils](https://github.com/indexofire/vpautils) — `OAgc.gbk` (32
O-antigen gene cluster loci) and `CPSgc.gbk` (163 capsular polysaccharide loci) —
into ONE merged `kind: cluster` database. It is a PURE built-in bundled database in
the `ecoli_dec`/`lm_doumith` pattern: no provider, no source URLs, no
`gapit db fetch vp_ok`. The committed three-scheme typing asset
(`src/gapit/data/dbs/vp_ok/typing.json`, generated reproducibly by
`scripts/generate_vp_typing.py`) installs with it, so screening renders
`gapit.cluster/2` with per-file phenotypes `o_group`, `k_group`, and the composed
`serotype` `"{o_group}:{k_group}"` (fallbacks `OUT`/`KUT`) — the O:K serotype from a
single scan:

```console
$ gapit screen RIMD.fna --db vp_ok --format json
{"schema": "gapit.cluster/2", ..., "phenotypes": {"o_group": {"phenotype": "O3"},
 "k_group": {"phenotype": "K6"}, "serotype": {"phenotype": "O3:K6"}}}
```

The database content is **MIT** (Bakta-annotated, via vpautils), so unlike the
kaptive providers it is license-clean for redistribution — a resolved snapshot ships
inside the wheel as one of the bundled databases (8.3 MiB committed artifacts, within
the 30 MiB budget), and first use materializes it with zero network. The
typing rules recompute sequence-unique key genes per locus (variants fold: OL4V1 →
O4, KL20V1 → K20; untyped OLU*/KLU* loci get no rule and fall back to OUT/KUT), and
the asset validates against the loci at build and materialization time — a drifted
upstream GBK fails `TYPING_UNKNOWN_GENE` instead of installing silently. Calibrated
against real genomes with RIMD 2210633 locked to O3:K6 (see the rule notes for the
documented per-locus tuning).

## Checking database freshness

`gapit db outdated` reports every installed database's age against the staleness
threshold. Against a fully installed datadir (output trimmed to three of twenty-three rows):

```console
$ gapit db outdated --days 30
NAME	FETCHED_AT	AGE_DAYS	STATUS
argannot	2026-09-17T23:22:06Z	2.58	ok
bacmet2	2026-09-17T23:13:51Z	2.58	ok
card	2026-09-17T23:14:25Z	2.58	ok
...
```

A database installed long ago simply reports `stale`:

```console
$ gapit db outdated --datadir /tmp/opencode/gapit-outdated-demo
NAME	FETCHED_AT	AGE_DAYS	STATUS
card	2020-01-01T00:00:00Z	2454.55	stale
```

| Status | Meaning |
|---|---|
| `ok` | Within the threshold |
| `stale` | `age_days` past `--days` (default 90; `--days 0` marks everything stale) |

Staleness is a report, never an error state: the command exits 0 however stale things are.
Exit 4 (`DATADIR_NOT_FOUND`, `DATADIR_EMPTY`) covers a missing datadir or one with no installed
databases; an unparseable `fetched_at` is `MANIFEST_MALFORMED` (exit 5). For bundled
databases the report is manifest-based as everywhere else: a materialized bundle ages
from its materialization time, and no bundled-vs-upstream comparison happens in v1
(see [Bundled databases](#bundled-databases-install-time-ready) — future work). For
agents, `--json` emits a `gapit.dboutdated/1` document (introspectable via
`gapit schema dboutdated`):

```console
$ gapit db outdated --days 3 --json
{
  "schema": "gapit.dboutdated/1",
  "databases": [
    {
      "db": "argannot",
      "fetched_at": "2026-09-17T23:22:06Z",
      "age_days": 2.58,
      "status": "ok",
      "upstream_version": ""
    },
    {
      "db": "card",
      "fetched_at": "2026-09-17T23:14:25Z",
      "age_days": 2.58,
      "status": "ok",
      "upstream_version": ""
    }
  ]
}
```

## Searching records across databases

`gapit db search TERM` looks up genes in the `records.jsonl` truth store of every installed
database — no BLAST, just a fast scan. Matching is case-insensitive substring by default;
`--exact` switches to full-field equality:

```console
$ gapit db search ctx-m --limit 3
argannot	(Bla)blaCTX-M-1	X92506:63-938		(Bla)blaCTX-M-1	876
argannot	(Bla)blaCTX-M-10	AF255298:1-873		(Bla)blaCTX-M-10	873
argannot	(Bla)blaCTX-M-100	FR682582:1-876		(Bla)blaCTX-M-100	876
$ gapit db search "blaCTX-M-1" --field gene --exact
ncbi	blaCTX-M-1	NG_048897.1	CEPHALOSPORIN	extended-spectrum class A beta-lactamase CTX-M-1	876
$ gapit db search virulence --db vfdb --field function --limit 2
vfdb	AAA92657	AAA92657	virulence	(AAA92657) unknown protein [TraJ (VF0241) - Invasion (VFC0083)] [Escherichia coli]	606
vfdb	AAC38364	AAC38364	virulence	(AAC38364) Orf1 [Ler (VF0189) - Regulation (VFC0301)] [Escherichia coli O127:H6 str. E2348/69]	390
```

Rows are `DB\tGENE\tACCESSION\tFUNCTION\tPRODUCT\tLENGTH`, streamed in database-then-file
order; `FUNCTION` joins the record's function classes with `;`. `--field` picks
`gene|accession|function|product|any` (default `any` searches all of them, one function class
at a time). `--limit N` caps the output (default 100; `0` = unlimited) and stderr notes a
truncation (`--quiet` silences it); zero hits exit 0 with empty stdout. `--json` prints one
JSON object per hit with the same fields as snake_case (`function` as an array):

```console
$ gapit db search "tet(M)" --field gene --exact --json | head -1
{"db":"card","gene":"tet(M)","accession":"AB039845.1:25-1945","function":["tetracycline"],"product":"Tet(M) is a ribosomal protection protein that confers tetracycline resistance. It is found on transposable DNA elements and its horizontal transfer between bacterial species has been documented.","length":1920}
```

An installed database whose `records.jsonl` is missing is skipped with a stderr warning during
a full scan, but targeting it explicitly (`--db NAME`) fails with exit 4 `DB_INCOMPLETE`; an
unknown `--db NAME` is a usage error (exit 2) listing what is installed.

## Providers

Nineteen providers ship with gapit. Every one of them — `card` and `vfdb` included —
downloads from its upstream source at fetch time: with the exception of the eight
audited, permissively licensed bundles above, nothing is bundled inside the package,
because several upstream licenses (CARD's McMaster terms, VFDB's CC BY-NC, Kaptive's
GPL-3.0) forbid redistribution inside an MIT-licensed distribution. The seven kaptive
providers are **cluster** databases (built through the cluster pipeline at
fetch time — see [Kaptive providers](#kaptive-providers-gpl-downloaded-on-fetch));
the only other cluster database, the typed `vp_ok`, is bundled-only
([vp_ok](#vp_ok-v-parahaemolyticus-ok-mit-typed-on-install)).

| Name | Maintainer | Content | dbtype |
|---|---|---|---|
| `ncbi` | NCBI | NCBI AMRFinderPlus (reference finder) curated AMR | nucl |
| `card` | McMaster University | CARD protein homolog resistance models | nucl |
| `resfinder` | DTU CGE | CGE ResFinder acquired resistance genes | nucl |
| `argannot` | IHU Méditerranée-Infection | ARG-ANNOT acquired resistance genes | nucl |
| `plasmidfinder` | DTU CGE | CGE PlasmidFinder replicons | nucl |
| `megares` | MEG Lab | MEGARes antimicrobial resistance genes | nucl |
| `ecoh` | Holt lab (srst2) | E. coli O and H antigens (srst2 EcOH) | nucl |
| `vfdb` | USTC (VFDB) | VFDB virulence factors (set A, nucleotide) | nucl |
| `ecoli_vf` | PHAC-NML | E. coli virulence factors (phac-nml) | nucl |
| `bacmet2` | University of Gothenburg | BacMet2 experimentally confirmed biocide/resistance genes (protein) | prot |
| `victors` | University of Chicago | Victors virulence factors | nucl |
| `upec_expec_vf` | FordeGenomics | UPEC/ExPEC virulence genes (FordeGenomics) | nucl |
| `kpsc_k` | Kaptive (klebgenomics) | K. pneumoniae species complex K locus (Kaptive; cluster, GPL-3.0, downloaded on fetch) | nucl |
| `kpsc_o` | Kaptive (klebgenomics) | K. pneumoniae species complex O locus (Kaptive; cluster, GPL-3.0, downloaded on fetch) | nucl |
| `kosc_k` | Kaptive (klebgenomics) | K. oxytoca species complex K locus (Kaptive; cluster, GPL-3.0, downloaded on fetch) | nucl |
| `kosc_o` | Kaptive (klebgenomics) | K. oxytoca species complex O locus (Kaptive; cluster, GPL-3.0, downloaded on fetch) | nucl |
| `ab_k` | Kaptive (Kenyon lab) | A. baumannii K locus (Kaptive; cluster, GPL-3.0, downloaded on fetch) | nucl |
| `ab_o` | Kaptive (Kenyon lab) | A. baumannii OC locus — official keyword `ab_o` (Kaptive; cluster, GPL-3.0, downloaded on fetch) | nucl |
| `ecoli_kps` | Kaptive (Gladstone lab) | E. coli group 2+3 capsular polysaccharide loci (Kaptive; cluster, GPL-3.0, downloaded on fetch) | nucl |

Protein databases (`bacmet2`) screen through `blastx`; nucleotide ones through `blastn`.
Cluster databases screen through the minimap2 cluster engine
([screen.md](./screen.md#cluster-databases-kind-cluster)).

### Kaptive providers (GPL, downloaded on fetch)

The seven kaptive providers wrap the reference databases of
[Kaptive](https://github.com/klebgenomics/Kaptive) (Wyres et al., J Clin Microbiol 2020 —
please cite Kaptive when you use results from these databases). The database NAMEs are the
**official install keywords** from the Kaptive v3 database docs
([Available databases](https://klebgenomics.github.io/Kaptive/db/overview.html#available-databases)),
and each fetches its raw GenBank file from the head (`main`) of the actively curated
per-species upstream repository — `klebgenomics/KpSC_surface_antigen_loci` (`kpsc_k`,
`kpsc_o`), `klebgenomics/KoSC-surface-antigen-loci` (`kosc_k`, `kosc_o`),
`johannajkenyon/Abaumannii_surface_polysaccharide_loci` (`ab_k`, `ab_o`), and
`rgladstone/EC-K-typing` (`ecoli_kps`) — so a fetched database tracks current upstream
curation rather than a frozen release. Those repos also ship `.toml` metadata with upstream
identity thresholds; gapit does not fetch it in v1 (a future `typing.json` source).

The database content is **GPL-3.0**, while gapit is MIT — so nothing kaptive is ever
bundled: `gapit db fetch kpsc_k` (and the other six keywords) downloads the GenBank file at
fetch time and builds a `kind: cluster` database through the cluster pipeline. The manifest
records the upstream URL, the `GPL-3.0 (database content)` license, and the citation note.
No typing model ships with them (phenotype calls stay null); you get kaptive-style best
locus calls, and you can later install your own `typing.json` semantics by rebuilding with
`gapit db build --typing`.

## Fetching databases

### The default set

`gapit db fetch all` installs the default set, `card` then `vfdb`. Like every
provider it downloads from upstream, transforms the records, and builds `sequences` plus
the BLAST index locally — so the index always matches the BLAST version actually
installed on your machine. Network access is required. A bare `gapit db fetch` prints the
command help: a multi-database download is always named explicitly (`all` or a single
NAME), never implied by an omitted argument.

### Fetching a named database

`gapit db fetch <name>` runs the full pipeline for one database. Here is `card`
(stderr progress lines, then a one-line JSON receipt on stdout):

```console
$ gapit db fetch --datadir /tmp/opencode/gapit-dbs-demo card
gapit: downloaded 1 source file(s)
gapit: read 6059 records from card
gapit: generated /tmp/opencode/gapit-dbs-demo/card/sequences
gapit: self-check passed for card
gapit: BLAST index built (nucl)
{"db":"card","records":6059,"dbtype":"nucl","destination":"/tmp/opencode/gapit-dbs-demo/card"}
```

| Receipt field | Type | Meaning |
|---|---|---|
| `db` | string | Database name (NAME) |
| `records` | integer | Records written to `records.jsonl` |
| `dbtype` | string | `nucl` or `prot` |
| `destination` | string | Installed database directory |

Every database downloads from its upstream source at fetch time, so fetches need network
access. Unknown names are rejected before the datadir is even resolved:

```console
$ gapit db fetch nosuchdb
{"schema":"gapit.error/1","code":"USAGE_ERROR","message":"unknown database: nosuchdb (available: argannot, bacmet2, card, ecoh, ecoli_vf, megares, ncbi, plasmidfinder, resfinder, upec_expec_vf, vfdb, victors)","context":{"db":"nosuchdb"}}
```

### Refetching

An existing database directory is never overwritten silently. Refetching without flags
fails with exit 4:

```console
$ gapit db fetch --datadir /tmp/opencode/gapit-dbs-demo card
{"schema":"gapit.error/1","code":"DB_ALREADY_EXISTS","message":"won't overwrite existing database card (use --force)","context":{"db":"card"}}
```

`--force` deletes and rebuilds in place:

```console
$ gapit db fetch --datadir /tmp/opencode/gapit-dbs-demo --force card
{"db":"card","records":6059,"dbtype":"nucl","destination":"/tmp/opencode/gapit-dbs-demo/card"}
```

## Installing a local file

`gapit db install SOURCE --sha256 HASH --output TARGET` is a verified local-file install.
It knows nothing about the provider catalog, archives, or sequence content: it streams SOURCE through
SHA256, compares the digest against `--sha256`, and only then atomically replaces TARGET.
A failed check leaves any existing TARGET untouched.

```console
$ sha256sum card.json
65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025  card.json
$ gapit db install card.json \
    --sha256 65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025 \
    --output /tmp/opencode/card-copy.json
{"destination":"/tmp/opencode/card-copy.json","sha256":"65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025"}
```

A wrong digest aborts with exit 5:

```console
$ gapit db install card.json \
    --sha256 0000000000000000000000000000000000000000000000000000000000000000 \
    --output /tmp/opencode/card-bad.json
{"schema":"gapit.error/1","code":"CHECKSUM_MISMATCH","message":"SHA256 mismatch for card.json","context":{"source":"card.json","expected":"0000000000000000000000000000000000000000000000000000000000000000","actual":"65838c8d4f160923fbf8c296c0f986bd5b55f7f5b5a15a7bbeb35412b945e025"}}
```

## What's inside a gapit-built database

```
<datadir>/<name>/
  records.jsonl         truth source: one Record JSON object per line
  sequences             generated FASTA projection (gapit/v1 tagged headers)
  sequences.n*|p*       BLAST index built from sequences
  gapit-manifest.json   provenance sidecar, written last
  typing.json           optional gapit.typing/1 or /2 scoring spec (typed dbs)
  floors.json           optional gapit.floors/1 per-gene identity floors (reads mode)
```

### records.jsonl

This is the editable truth. Every screening artifact is generated from it, so editing
records and refetching with `--force` regenerates everything downstream. One JSON object
per line; the first card record, verbatim:

```json
{"db":"card","gene":"23S_rRNA_(adenine(2058)-N(6))-methyltransferase_Erm(A)","sequence":"ATGAAACAGAAAAACCCGAAAAATACGCAAAATTTCATTACATCTAAAAAGCATGTAAAGGAAATATTAAAATATACGAATATCAATAAACAAGATAAAATAATAGAAATTGGGTCAGGAAAAGGACATTTTACCAAGGAACTTGTGGAAATGAGTCAACGGGTGAATGCTATAGAGATTGATGAAGGTTTATGTCATGCCACGAAAAAAGCAGTTGAACCTTTTCAGAATATAAAAGTTATTCATGAGGATATTTTGAAGTTTAGCTTTCCTAAAAATACAGACTATAAAATATTTGGTAATATTCCCTACAATATTAGTACTGATATTGTAAAAAAGATTGCTTTTGATAGTCAAGCGAAATATAGCTACCTTATTGTAGAGAGGGGATTTGCTAAAAGGTTGCAAAATACCCAACGAGCTTTAGGTTTGCTGTTAATGGTGGAAATGGATATAAAAATTCTTAAAAAAGTGCCACGAGCATATTTTCACCCTAAGCCTAATGTAGATTCTGTATTGATTGTACTTGAAAGGCATAAACCATTTATTTTAAAGAAGGACTACAAAAAGTATAGATTTTTCGTTTATAAATGGGTAAACAGGGAATATCATGTTCTTTTTACTAAAAATCAATTAAGACAGGTGCTGAAGCATGCGAATGTTACTGATCTTGATAAATTATCCAATGAACAATTTTTGTCTGTTTTCAATAGTTACAAATTATTTCAATAA","accession":"AF002716.1:210-942","function":["lincosamide","macrolide","streptogramin"],"product":"Variant of ErmA (ARO:3000347) found in Streptococcus pyogenes. Confers the MLSb phenotype.","source_id":"3005099"}
```

| Field | Type | Meaning |
|---|---|---|
| `db` | string | Source database name |
| `gene` | string | Gene name, the identity gapit reports on |
| `sequence` | string | Normalized sequence (uppercase; `N`/`X` for ambiguous bases) |
| `accession` | string | Accession or coordinate span, empty when unpublished |
| `function` | array of string | Functional categories, sorted (see the vocabulary below) |
| `product` | string | Product description, `n/a` when absent |
| `source_id` | string | Provider-native record identifier |

### gapit-manifest.json

The manifest certifies the build: what was fetched, when, from where, and with which tool
versions. A real one, from the plasmidfinder database:

```json
{
  "schema": "gapit.manifest/1",
  "name": "plasmidfinder",
  "source_urls": [
    "https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/HEAD.zip"
  ],
  "fetched_at": "2026-09-17T23:11:36Z",
  "sha256": "a18dc9dab762fe1e23a90c11314ee179c37c59dad833a5f465a78dc245970c2a",
  "n_records": 488,
  "dbtype": "nucl",
  "header_format": "gapit/v1",
  "upstream_version": "",
  "tool": {
    "name": "gapit",
    "version": "0.6.4"
  },
  "makeblastdb_version": "blastn: 2.17.0+",
  "minimap2_version": "2.31-r1302"
}
```

| Field | Type | Meaning |
|---|---|---|
| `schema` | string | Always `gapit.manifest/1` |
| `name` | string | Database name |
| `source_urls` | array of string | Upstream URLs the content came from |
| `fetched_at` | string | ISO-8601 UTC timestamp, second precision |
| `sha256` | string | Digest of the built `sequences` file, 64 lowercase hex chars |
| `n_records` | integer | Record count, matches `records.jsonl` |
| `dbtype` | string | `nucl` or `prot`, chosen explicitly at build time |
| `header_format` | string | Always `gapit/v1` for gapit-built databases |
| `upstream_version` | string | Upstream release label when the source has one |
| `license` | string | Database-content license, present only when the provider declares one (kaptive: `GPL-3.0 (database content)`; card: McMaster non-commercial terms; vfdb: `CC BY-NC 4.0`) |
| `note` | string | Free-text provenance note (kaptive: source file + citation), present only when set |
| `tool` | object | `{name, version}` of the gapit that built the database |
| `makeblastdb_version` | string | BLAST+ version that built the index |
| `minimap2_version` | string | minimap2 version in the build environment (reads mode indexes in memory; no `.mmi` is built) |

`records.jsonl` and the manifest are file contracts. They never appear on stdout and are
not registered with `gapit schema` (the stdio/listing documents `dblist`, `dboutdated`,
and `datadir` are).

### floors.json — per-gene identity floors

An optional sidecar declaring a MINIMUM alignment identity (%) per gene for reads-mode
presence (schema `gapit.floors/1`, introspect with `gapit schema floors`; see
[reads.md](./reads.md#per-gene-identity-floors-gapitfloors1-database-side) for the
screening semantics and the `pic` worked example):

```json
{"schema": "gapit.floors/1", "default": null, "genes": {"pic": 90.0}}
```

| Field | Type | Meaning |
|---|---|---|
| `schema` | string | Always `gapit.floors/1` |
| `default` | number or null | Floor for genes not listed in `genes` (`null` = no floor) |
| `genes` | object | `{gene: minimum identity %}`, values in [0, 100] |

The file's presence IS the flag — no manifest field, no CLI switch. A database without
the sidecar screens byte-identically to before (the floors are a db-driven opt-in, not a
CLI contract change), and the blastn contig path never reads them. Built-in via
`db build --floors` (next sections), or drop the file into an existing db directory; it
is validated at screening time too (malformed content fails with `FLOORS_MALFORMED`,
exit 4). The bundled `ecoli_dec` ships `{"genes": {"pic": 90.0}}` — the SPATE-homolog
false-positive fix.

## Header formats and compatibility

gapit reads two header formats and detects them per record:

- **abricate legacy**: `>DB~~~GENE~~~ACCESSION~~~RESISTANCE PRODUCT`. Fields may be
  missing; missing trailing fields decode to empty strings.
- **gapit-native `gapit/v1`**: `>gapit|db=<v>|gene=<v>|acc=<v>|func=<v> PRODUCT` with
  percent-encoded values, so gene names and products survive BLAST and minimap2 byte for
  byte, spaces and pipes included.

Both decode into the same four identity fields, and the TSV/JSON outputs are identical in
shape either way. For native databases the `RESISTANCE` column carries the `func` value.

> **Warning:** abricate cannot read gapit-built databases. The `gapit/v1` header format is
> not the `~~~` convention. gapit reading abricate datadirs is one-directional by design;
> keep a legacy copy if a downstream tool still runs abricate.

### The `func` vocabulary

Each native provider fills `function` from a fixed vocabulary, and these values flow
through to the TSV `RESISTANCE` column and the JSON `resistance` field:

| Provider(s) | `func` value(s) |
|---|---|
| `ncbi`, `resfinder`, `argannot`, `card` | The source's antibiotic classes |
| `megares` | The MEGARes class field |
| `ecoh` | `H-antigen`, `O-antigen`, or an antigen fallback derived from the allele |
| `vfdb`, `ecoli_vf`, `victors`, `upec_expec_vf` | `virulence` |
| `plasmidfinder` | `replicon` |
| `bacmet2` | `biocide` |

## Bringing your own database

`gapit db build NAME FASTA` turns any FASTA of reference genes into a fully built
gapit-native database — `records.jsonl`, the `sequences` projection with `gapit/v1`
headers, the BLAST index, and the manifest — in one command.
Here is a two-gene synthetic FASTA plus a metadata TSV (fields below), run against a
scratch datadir:

```console
$ gapit db build tinyamr my_genes.fa --datadir ./db --tsv my_meta.tsv
gapit: generated /tmp/opencode/gapit-build-demo/db/tinyamr/sequences
gapit: self-check passed for tinyamr
gapit: BLAST index built (nucl)
{"db":"tinyamr","records":2,"dbtype":"nucl","destination":"/tmp/opencode/gapit-build-demo/db/tinyamr"}
```

stderr carries the per-step progress, stdout the one-line JSON receipt (same fields
as `db fetch`, see the table above). The database screens immediately:

```console
$ gapit screen --datadir ./db --db tinyamr contig.fa
Processing: contig.fa
Found 1 genes in contig.fa
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
contig.fa	contig1	1	240	+	syn_betalac	1-240/240	===============	0/0	100.00	100.00	tinyamr	SYN-0001	synthetic class A beta-lactamase	ampicillin;cephalosporin
```

The `ACCESSION` and `RESISTANCE` values came from the TSV merge; the `PRODUCT` from
the FASTA description. `--dbtype nucl|prot` forces the molecule type; by default the
abricate mol-type heuristic decides from the sequences themselves. Input may be plain,
`.gz`, or `.bz2`.

### Header detection

The header kind is detected per record, so mixed files work:

| FASTA header | Gene | Accession | Function | Product |
|---|---|---|---|---|
| `>syn_betalac synthetic class A beta-lactamase` | `syn_betalac` | — | — | `synthetic class A beta-lactamase` |
| `>olddb~~~sul1~~~U12338.4:1-940~~~SULFONAMIDE sulfonamide resistance` | `sul1` | `U12338.4:1-940` | `SULFONAMIDE` | `sulfonamide resistance` |
| `>gapit\|db=old\|gene=sul2\|acc=X\|func=streptomycin aminoglycoside` | `sul2` | `X` | `streptomycin` | `aminoglycoside` |

Plain headers without any description text fall back to `--description TEXT`, then
to the gene name. The `db` field of every record is always NAME (the database being
built), and `source_id` keeps the original id token verbatim. Malformed `gapit|`
headers fail the build with `HEADER_MALFORMED` (exit 4); structurally invalid FASTA
fails with `INVALID_FASTA` (exit 5).

### Merging metadata from a TSV

`--tsv FILE` adds or overrides accession and function classes per gene. The rules:

- The header row is mandatory and must contain a `gene` column; `accession` and
  `function` are optional per file (absent columns are simply not merged), and extra
  columns are ignored. A missing `gene` column fails with exit 5
  `METADATA_MALFORMED`.
- Rows are keyed by gene: the **first row wins** on duplicates, with a warning on
  stderr (`--quiet` silences it).
- The `function` column is `;`-separated for multiple classes, e.g.
  `ampicillin;cephalosporin`.
- Genes present only in the TSV produce a stderr warning and are skipped — the FASTA
  is the truth for what exists.

The example TSV used above, in full:

```console
$ cat my_meta.tsv
gene	accession	function
syn_betalac	SYN-0001	ampicillin;cephalosporin
```

### Per-gene identity floors (`--floors`)

`--floors FILE` installs a `gapit.floors/1` document as the database's `floors.json`
sidecar — the per-gene minimum alignment identity for reads-mode presence (see
[reads.md](./reads.md#per-gene-identity-floors-gapitfloors1-database-side)). The
document is validated before any artifact is written: every gene in `genes` must exist
in the FASTA (`FLOORS_UNKNOWN_GENE`, exit 4 — a floor for a missing gene is silently
dead safety config), values must lie in [0, 100] and the structure parse
(`FLOORS_MALFORMED`, exit 4). The file is then copied in byte-identical; the file's
presence is the flag (no manifest field changes). Gene builds only — a GBK/GFF cluster
input rejects `--floors` as a usage error.

```console
$ cat my_floors.json
{"schema": "gapit.floors/1", "default": null, "genes": {"syn_betalac": 95.0}}
$ gapit db build tinyamr my_genes.fa --datadir ./db --tsv my_meta.tsv --floors my_floors.json
{"db":"tinyamr","records":2,"dbtype":"nucl","destination":"/tmp/opencode/gapit-build-demo/db/tinyamr"}
$ cat ./db/tinyamr/floors.json
{"schema": "gapit.floors/1", "default": null, "genes": {"syn_betalac": 95.0}}
```

### Rebuilding

Like `db fetch`, an existing database is never overwritten silently:

```console
$ gapit db build tinyamr my_genes.fa --datadir ./db
{"schema":"gapit.error/1","code":"DB_ALREADY_EXISTS","message":"won't overwrite existing database tinyamr (use --force)","context":{"db":"tinyamr"}}
$ gapit db build tinyamr my_genes.fa --datadir ./db --tsv my_meta.tsv --force
gapit: generated /tmp/opencode/gapit-build-demo/db/tinyamr/sequences
gapit: self-check passed for tinyamr
gapit: BLAST index built (nucl)
{"db":"tinyamr","records":2,"dbtype":"nucl","destination":"/tmp/opencode/gapit-build-demo/db/tinyamr"}
```

Because `records.jsonl` is the editable truth, editing it (or the FASTA) and
rebuilding with `--force` regenerates every downstream artifact.

### Renaming an installed database

`gapit db rename OLD NEW` re-homes a database under a new name without a re-download:

```console
$ gapit db rename cam cam_vf
gapit: renaming cam -> cam_vf (gene database)
gapit: BLAST index built (nucl)
{"db":"cam_vf","renamed_from":"cam","records":17,"kind":"gene","destination":".../cam_vf"}
```

Gene databases carry the name inside every `records.jsonl` record and every `gapit/v1`
header, so the directory is rebuilt through the standard build pipeline with `db`
re-pointed (records, `sequences`, self-check, BLAST index, manifest). Cluster databases
carry the name only in the manifest — locus headers are bare ids — so their content is
copied verbatim and the manifest re-stamped. Typing/floors/feature sidecars travel along,
and both paths are staged: a failure leaves OLD untouched, and OLD is removed only after
NEW is published (the manifest provenance, including `fetched_at`, is preserved).

Use it when a name predates the [naming convention](#database-naming); it is also the
generic escape hatch for a lab's own naming choices. `--force` replaces an existing
target. Directories without a manifest (abricate-style imports) have no name to rewrite
and are refused with `DATABASE_NOT_FOUND` — rebuild them with `db build` instead.

### The manual (abricate-style) path

`db build` is the recommended route, but any abricate-format datadir also works
as-is: point `--datadir` at a directory containing `<name>/sequences` (a nucleotide
FASTA with `~~~` headers) and run `gapit setupdb` once to build the BLAST index. The
mol-type heuristic from abricate picks `nucl` vs `prot`.

End to end with a tiny three-gene database copied into a scratch datadir:

```console
$ ls /tmp/opencode/gapit-own-db/db/tinyamr
sequences
$ gapit screen --datadir /tmp/opencode/gapit-own-db/db --db tinyamr /tmp/opencode/gapit-own-db/gap.fa
Processing: /tmp/opencode/gapit-own-db/gap.fa
{"schema":"gapit.error/1","code":"DATABASE_NOT_INDEXED","message":"Database /tmp/opencode/gapit-own-db/db/tinyamr/sequences is not indexed, please try: gapit setupdb","context":{"db":"/tmp/opencode/gapit-own-db/db/tinyamr/sequences"}}
$ gapit setupdb --datadir /tmp/opencode/gapit-own-db/db
Indexed tinyamr (3 sequences, nucl)
$ gapit screen --datadir /tmp/opencode/gapit-own-db/db --db tinyamr /tmp/opencode/gapit-own-db/gap.fa
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
/tmp/opencode/gapit-own-db/gap.fa	contig1	1	97	+	sul1	1-94/94	========/======	1/3	100.00	96.91	tinyamr	U12338.4:1-940	sulfonamide-resistant dihydropteroate synthase Sul1	SULFONAMIDE
```

`gapit setupdb` indexes every subdirectory with a readable `sequences` file and exits 0
when all of them have `sequences.nin` or `sequences.pin`.

### Reads mode and legacy datadirs

For reads screening (see [reads.md](./reads.md)), minimap2 always indexes the `sequences`
FASTA in memory; no `.mmi` is built, and one left in a datadir by an older gapit is never
consulted. Persisted indexes were tried and rejected: a default-built `.mmi` overrides the
`-x sr` preset's indexing parameters (minimap2 warns "-k, -w or -H overridden by prebuilt
index"), which misassigns close homologs (CTX-M/SHV allele divergence in the deciding
benchmark) and ran slower than in-memory indexing at current database scale. Legacy
abricate datadirs take the same in-memory path.

A complete `db build` walkthrough with worked examples for every header format, protein
databases, metadata TSVs, and a troubleshooting table lives in
[Custom databases](./custom-db.md).

## Cluster databases (GBK/GFF)

`gapit db build NAME input.gbk` (or `.gbff`, `.gb`, `.gff`, `.gff3`, each plain or
`.gz`/`.bz2`) builds a different kind of database: one record = one gene locus. The input
kind is detected by suffix (`--kind gene|cluster` must agree or the build is a usage
error); the FASTA gene pipeline of the previous section is untouched.

```console
$ gapit db build mycps cps_loci.gbk --datadir ./db
gapit: generated ./db/mycps/sequences
gapit: self-check passed for mycps
gapit: BLAST index built (nucl)
{"db":"mycps","records":2,"dbtype":"nucl","destination":"./db/mycps"}
```

Parsing (no Biopython): GenBank records become loci (`LOCUS` name = locus id) with their
`CDS` features as genes; the source-feature notes carry the label and type (kaptive
`K locus:`/`K type:` style, with Bakta-style fallbacks to the locus id). GFF3 inputs take
the sequences from an embedded `##FASTA` block or a `<stem>.fa/.fna/.fasta` sidecar.
Compound `join()` CDS locations are rejected with a typed error — split the record or use
simple CDS annotations.

The artifacts differ from a gene database: `sequences` holds one bare-locus-id record per
locus (plain headers — abricate cannot read these), `features.json` is the
`gapit.features/1` feature table (introspect with `gapit schema features`), the manifest
declares `kind: cluster`, and there is no `records.jsonl`. Screening dispatches to the
minimap2 cluster engine ([screen.md](./screen.md#cluster-databases-kind-cluster)).

### typing.json — declarative phenotype scoring

`--typing FILE` installs a validated `gapit.typing/1` or `gapit.typing/2` document into
the database as `typing.json`; screening then annotates every best call with a phenotype
and an explainable score breakdown (the PHENOTYPE TSV column and `phenotype_detail` in
gapit.cluster/1 —
[outputs.md](./outputs.md#gapitcluster1-cluster-database-screening)). A minimal
`weighted_genes` example:

```json
{
  "schema": "gapit.typing/1",
  "rules": [
    {
      "model": "weighted_genes",
      "phenotype": "K1",
      "weights": {"wzx": 1.0, "wzy": 3.0},
      "negative": {"rfaD": -1.0},
      "identity_floor": 95.0,
      "require_any": ["wzx", "wzy"]
    }
  ],
  "cutoff": 0.9,
  "ambiguity_margin": 0.05,
  "fallback": "unknown"
}
```

Three rule kinds exist — `weighted_genes` (gene presence with an identity floor, negative
markers, and an any-of gate), `cluster_match` (weighted locus coverage/identity/key-genes
components with per-component floors), and `learned_linear` (a trained sigmoid over named
features `gene:<id>:<cov|ident|present>` / `cluster:<locus>:<coverage|identity>`) — plus
the document-level `cutoff`, `ambiguity_margin`, and `fallback` the decision layer applies.
The full schema prints with `gapit schema typing`. Validation is strict at build time: a
malformed document fails with `TYPING_MALFORMED`, and a rule referencing a gene or locus
the input does not carry fails with `TYPING_UNKNOWN_GENE` before any artifact is written
(the same check reruns at screen time, so hand-edited databases cannot smuggle dead
references).

`gapit.typing/2` wraps the same rules in NAMED schemes (`schemes:
[{"name": "pathotype", "rules": [...], "cutoff": ..., "ambiguity_margin": ...,
"fallback": ...}]`, names `[a-z0-9_]` and unique); a `/1` document loads as one
anonymous `default` scheme, so existing typed databases keep their output. `--typing`
also works on GENE (FASTA) builds — `weighted_genes` rules only in v1, gene references
checked against the FASTA records. A typed gene db screens exactly like an untyped one
(pure gene detection); its designations come from the two-stage pipeline — `gapit screen
-o result.tsv` writes the table, and `gapit typing result.tsv` renders one call per
scheme as `gapit.typing_result/1`
([typing.md](./typing.md#the-two-stage-designation-workflow)).

typing/2 stage 2 adds six rule/scheme primitives (all optional and additive; full
detail in `gapit schema typing`): the `exact_set` rule (Doumith/Shigella-style marker
tables — score 1.0 iff every `requires` gene is present AND every `excludes` gene is
absent, floors defaulting to 90/90; a 1.0 tie involving a satisfied exact_set resolves
by declaration order instead of the ambiguity margin), an optional `coverage_floor` on
`weighted_genes`/`exact_set`, scheme-level `control_gene` (prs/ipaH-style gate: when
absent the whole scheme outputs its fallback with a "control gene absent" note),
`unique_group` + `mixed_phenotype` (when more than one gene of a group is present, the
scheme calls the mixed phenotype with the pair in `ambiguous`), `compose` schemes
(`"{o_group}:{k_group}"` rendered from sibling schemes' calls — fallback strings flow
through, an ambiguous ingredient yields a null composition carrying that ingredient's
variants), and rule-level `notes` surfaced verbatim on winning calls.

Six researched designation schemes — the complete Doumith Listeria table, ShigaTyper-semantics
Shigella/EIEC, the meningotype serogroup panel with allele probes, the Vibrio parahaemolyticus
O/K Kaptive pattern, V. cholerae O1/O139 + Ogawa/Inaba, and a documented DEC placeholder — are
encoded as validated example documents (synthetic marker fixtures, self-checked by the test
suite) and walked through in [Typing schemes](./typing.md).

Tuning a typing document against labeled assemblies is what
`scripts/cluster_calibration.py` is for (developer tool): it screens every labeled sample,
prints the per-expected-phenotype score distribution, the called×expected agreement
matrix, and the divergence list — the loop a future `learned_linear` training flow will
consume:

```console
$ python scripts/cluster_calibration.py mycps --datadir ./db labels.tsv
```
