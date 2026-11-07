# Quickstart

A complete first session, offline and reproducible, using the fixture databases that
ship with the repo and wheel. You screen two contig files, read the results as TSV and
JSON, type a sample against the bundled DEC panel, and summarize reports into one
matrix.

Run everything from the repo root with the gspec environment on PATH
(`export PATH="$PWD/.pixi/envs/default/bin:$PATH"`, or prefix every command with
`pixi run`).

## 1. Set up a throwaway datadir

The fixture db at `tests/data/db/tinyamr` holds three short AMR genes: `tetA`,
`blaTEM-1`, and `sul1`. Copy it into a temp datadir and let `gspec db setup` bootstrap
the rest: it materializes the seven bundled databases that ship in the wheel (zero
network — including `eco_dec`, which you will type with in section 4) and indexes
everything, no `makeblastdb` by hand:

```bash
mkdir -p /tmp/gspec-quickstart/db
cp -r tests/data/db/tinyamr /tmp/gspec-quickstart/db/
export GSPEC_DATADIR=/tmp/gspec-quickstart/db
gspec db setup
```

All progress goes to stderr; one run prints (trimmed):

```console
gspec: materializing bundled database cje_cps (42 records) into /tmp/gspec-quickstart/db
gspec: materializing bundled database eco_dec (16 records) into /tmp/gspec-quickstart/db
gspec: materializing bundled database lmo_sero (5 records) into /tmp/gspec-quickstart/db
gspec: materializing bundled database vpa_ok (195 loci) into /tmp/gspec-quickstart/db
...
gspec: indexed eco_dec (16 sequences, nucl)
gspec: indexed tinyamr (3 sequences, nucl)
```

`gspec db list` shows the database catalog and install state (trimmed below; `--json`
returns the `gspec.dblist/1` document). The fixture `tinyamr` is a plain custom database,
not a catalog database, so it does not appear there — the screen in the next section confirms it
is usable:

```console
$ gspec db list
NAME	DESCRIPTION	STATUS	DBTYPE	KIND	TYPING	PROVIDER
argannot	ARG-ANNOT acquired resistance genes	available	nucl	gene	-	IHU Méditerranée-Infection
cje_cps	Campylobacter jejuni Penner capsule (CPS) multiplex PCR panel	installed (42)	nucl	gene	command	gspec-curated (Poly 2015; CC0 + NCBI public domain)
```

(twenty-eight catalog databases in total: the seven bundled ones now read `installed` in
this datadir, the rest read `available` until you `gspec db fetch` them)

## 2. Screen a contig file

```console
$ gspec screen tests/data/contigs/full.fa --db tinyamr
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
tests/data/contigs/full.fa	contig1	1	79	+	tetA	1-79/79	===============	0/0	100.00	100.00	tinyamr	NC_000913.3:100-900	tetracycline efflux pump TetA	TETRACYCLINE
```

The TSV is abricate-compatible, drop it into any pipeline that already parses abricate.
`Processing:` progress lines go to stderr; stdout carries data only. The columns:

| Column | Meaning |
|---|---|
| `#FILE` | Input file, path as given on the command line |
| `SEQUENCE` | Contig holding the hit |
| `START`, `END` | Hit coordinates on the contig (1-based) |
| `STRAND` | Gene orientation, `+` or `-` |
| `GENE` | Gene name from the database header |
| `COVERAGE` | Alignment span over the gene, `range/length` |
| `COVERAGE_MAP` | Fixed-width coverage sketch of the gene; `=` covered, `/` marks a gap run |
| `GAPS` | Gap openings/gap bases in the alignment |
| `%COVERAGE` | Fraction of the gene covered (unrounded internally, shown to 2 decimals) |
| `%IDENTITY` | BLAST percent identity, shown to 2 decimals |
| `DATABASE` | Database screened against |
| `ACCESSION` | Upstream accession from the database header |
| `PRODUCT` | Gene product description |
| `RESISTANCE` | Resistance class(es) from the database header |

Defaults follow abricate: a hit survives with `--minid 80` and `--mincov 80`
(percent identity, percent gene coverage). See [Screening](./screen.md) for every flag
and the filtering rules.

## 3. The same result as JSON

```bash
gspec screen tests/data/contigs/full.fa --db tinyamr --format json
```

```json
{
  "schema": "gspec.report/1",
  "tool": {
    "name": "gspec",
    "version": "0.10.0"
  },
  "created_at": "2026-09-19T01:11:06Z",
  "params": {
    "db": "tinyamr",
    "minid": 80.0,
    "mincov": 80.0,
    "threads": 1
  },
  "files": [
    {
      "file": "tests/data/contigs/full.fa",
      "hits": [
        {
          "sequence": "contig1",
          "start": 1,
          "end": 79,
          "strand": "+",
          "gene": "tetA",
          "coverage": "1-79/79",
          "coverage_map": "===============",
          "gaps": "0/0",
          "coverage_pct": 100.0,
          "identity_pct": 100.0,
          "database": "tinyamr",
          "accession": "NC_000913.3:100-900",
          "product": "tetracycline efflux pump TetA",
          "resistance": "TETRACYCLINE"
        }
      ]
    }
  ]
}
```

The document carries `"schema": "gspec.report/1"`, so consumers can pin the contract;
`gspec schema report` prints its JSON Schema. Field-level detail lives in
[Outputs](./outputs.md).

## 4. Type your samples

Screening only detects genes. Designation (pathotype, serogroup) is a second command:
`gspec screen` writes the gene table, `gspec typing` reads it back and scores the
database's typing schemes. The bundled `eco_dec` panel carries the dual-scheme DEC
(diarrheagenic E. coli) designation. Use a test-fixture assembly as a stand-in for
yours:

```bash
cp tests/data/typing/dec_s3_stx2a_escV_aggR_uidA.fasta /tmp/gspec-quickstart/sample.fna
gspec screen /tmp/gspec-quickstart/sample.fna -d eco_dec -o /tmp/gspec-quickstart/result.tsv -p -q
gspec typing /tmp/gspec-quickstart/result.tsv -q
```

```console
FILE	SCHEME	PHENOTYPE	GENES	CONFIDENCE	SCORE	RUNNER_UP	NOTES
sample.fna	gb4789_6	EHEC	aggR;escV;stx2a;uidA	high	1.0000	EAEC (1.0000)	...
sample.fna	risk_monitoring	EHEC	aggR;escV;stx2a;uidA	high	1.0000	EAEC (1.0000)	...
```

One row per scheme: this isolate carries `aggR`, `escV`, `stx2a`, and `uidA`, so both
schemes call **EHEC** at full score (the `NOTES` column, trimmed above, spells out each
scheme's rule). The two stages also pipe straight into each other — screen's stdout is
typing's stdin:

```bash
gspec screen /tmp/gspec-quickstart/sample.fna -d eco_dec -p -q | gspec typing -q
```

Reads work the same way through the wildcard path: when every positional file is FASTQ,
samples are auto-paired from the filenames (`_R1`/`_R2`, `_1`/`_2`) and `--jobs`
parallelizes samples:

```bash
mkdir -p /tmp/gspec-quickstart/reads
cp tests/data/reads/s1_1.fq.gz tests/data/reads/s1_2.fq.gz \
  tests/data/reads/s2_R1.fq.gz tests/data/reads/s2_R2.fq.gz /tmp/gspec-quickstart/reads/
cd /tmp/gspec-quickstart/reads
gspec screen -d eco_dec *.fq.gz -j 4
```

```console
Screening sample s1 reads: s1_1.fq.gz, s1_2.fq.gz
Screening sample s2 reads: s2_R1.fq.gz, s2_R2.fq.gz
Detected 0 present genes in sample s2
Detected 0 present genes in sample s1
#SAMPLE	GENE	BREADTH%	DEPTH	READS	PRESENT	DATABASE	ACCESSION	PRODUCT
```

The demo reads carry tetX, not DEC markers, so the table stays header-only — zero hits
are normal output, not an error. Full designation semantics live in
[Typing schemes](./typing.md); the reads surface in [Screening reads](./reads.md).

## 5. Summarize reports into a matrix

Save two reports, then fold them into one gene-by-file matrix. The second fixture file,
`gap.fa`, carries `sul1` with gaps in the alignment:

```bash
gspec screen tests/data/contigs/full.fa --db tinyamr > full.tsv
gspec screen tests/data/contigs/gap.fa --db tinyamr > gap.tsv
gspec summary full.tsv gap.tsv
```

```console
#FILE	NUM_FOUND	sul1	tetA
full.tsv	1	-	+
gap.tsv	1	+	-
```

Each row is one report file; columns are the union of genes found. Cells hold the
presence call: `+` when the report has a hit for the gene, `.` when absent (pass
`--identity` or `--coverage` for the numbers). `NUM_FOUND` counts distinct genes. Add
`--format json|md` for machine- or human-readable matrices; details in [Summary](./summary.md).

## 6. Try the cluster engine

Gene **cluster** databases are a second database kind: instead of individual genes, one
record is a whole locus (Kaptive-style antigen loci, capsule clusters), built from
GenBank/GFF input. The fixture `tests/data/cluster/screening.gbk` holds two synthetic
loci; `--typing` installs a `gspec.typing/1` phenotype scoring spec alongside:

```bash
gspec db build tinykps tests/data/cluster/screening.gbk \
  --datadir /tmp/gspec-quickstart/db \
  --typing tests/data/cluster/typing_screen.json
# pull one locus out of the db as a one-locus query assembly
awk '/^>locusA$/{p=1} /^>/{if($0!~/>locusA$/)p=0} p' \
  /tmp/gspec-quickstart/db/tinykps/sequences > /tmp/gspec-quickstart/locusA.fa
gspec screen /tmp/gspec-quickstart/locusA.fa --db tinykps
```

Screening a cluster database dispatches to the minimap2 cluster engine; the TSV shape
changes to one best-locus call per file, and the typing spec turns the call into a
phenotype:

```console
Processing: /tmp/gspec-quickstart/locusA.fa
Best locus in /tmp/gspec-quickstart/locusA.fa: locusA
FILE	BEST_LOCUS	TYPE	PHENOTYPE	COVERAGE	IDENTITY	PRESENT	PARTIAL	MISSING_IDS
/tmp/gspec-quickstart/locusA.fa	locusA	KL101	K101	100.00	100.00	3	0	-
```

`PHENOTYPE` is `-` when no phenotype clears the spec's cutoff or two tie inside its
ambiguity margin (screen both loci at once and the K101/K102 tie does exactly that).
The seven kaptive providers (`gspec db fetch -d kpsc_k`, `kpsc_o`, `kosc_k`, `kosc_o`,
`ab_k`, `ab_o`, `ecoli_kps`) install
real Kaptive locus databases the same way. Full detail in
[Screening](./screen.md#cluster-databases-kind-cluster) and
[Databases](./databases.md#cluster-databases-gbkgff).

## Next steps

- [Screening](./screen.md): all contig-mode flags, thresholds, multiple inputs, `--jobs`
- [Screening reads](./reads.md): FASTQ input through minimap2
- [Typing schemes](./typing.md): the two-stage designation pipeline and the typing document format
- [Databases](./databases.md): install real databases (`gspec db fetch`), providers, datadirs, cluster databases and typing specs
- [Outputs](./outputs.md): formats, schemas, error envelopes, exit codes
- [MCP server](./mcp.md) and [Agent guide](./agents.md): driving gspec from agents
