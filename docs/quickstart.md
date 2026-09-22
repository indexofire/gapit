# Quickstart

A complete first session, offline and reproducible, using the tiny fixture database from
the test suite. You screen two contig files, read the results as TSV and JSON, and
summarize both reports into one matrix.

Run everything from the repo root with the gapit environment on PATH
(`export PATH="$PWD/.pixi/envs/default/bin:$PATH"`, or prefix every command with
`pixi run`).

## 1. Set up a throwaway datadir

The fixture db at `tests/data/db/tinyamr` holds three short AMR genes: `tetA`,
`blaTEM-1`, and `sul1`. Copy it into a temp datadir and build the BLAST index, exactly
what the test suite does for its MCP fixture:

```bash
mkdir -p /tmp/gapit-quickstart/db
cp -r tests/data/db/tinyamr /tmp/gapit-quickstart/db/
makeblastdb -in /tmp/gapit-quickstart/db/tinyamr/sequences \
  -title tinyamr -dbtype nucl -logfile /dev/null
export GAPIT_DATADIR=/tmp/gapit-quickstart/db
```

`gapit list` confirms the database is visible:

```console
$ gapit list
DATABASE	SEQUENCES	DBTYPE	DATE
tinyamr	3	nucl	2026-Sep-19
```

## 2. Screen a contig file

```console
$ gapit screen tests/data/contigs/full.fa --db tinyamr
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
gapit screen tests/data/contigs/full.fa --db tinyamr --format json
```

```json
{
  "schema": "gapit.report/1",
  "tool": {
    "name": "gapit",
    "version": "0.2.0"
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

The document carries `"schema": "gapit.report/1"`, so consumers can pin the contract;
`gapit schema report` prints its JSON Schema. Field-level detail lives in
[Outputs](./outputs.md).

## 4. Summarize reports into a matrix

Save two reports, then fold them into one gene-by-file matrix. The second fixture file,
`gap.fa`, carries `sul1` with gaps in the alignment:

```bash
gapit screen tests/data/contigs/full.fa --db tinyamr > full.tsv
gapit screen tests/data/contigs/gap.fa --db tinyamr > gap.tsv
gapit summary full.tsv gap.tsv
```

```console
#FILE	NUM_FOUND	sul1	tetA
full.tsv	1	.	100.00
gap.tsv	1	100.00	.
```

Each row is one report file; columns are the union of genes found. Cells hold
%COVERAGE (`.` when absent), `NUM_FOUND` counts distinct genes. Add `--format json|md`
for machine- or human-readable matrices; details in [Summary](./summary.md).

## Next steps

- [Screening](./screen.md): all contig-mode flags, thresholds, multiple inputs, `--jobs`
- [Screening reads](./reads.md): FASTQ input through minimap2
- [Databases](./databases.md): install real databases (`gapit db fetch`), providers, datadirs
- [Outputs](./outputs.md): formats, schemas, error envelopes, exit codes
- [MCP server](./mcp.md) and [Agent guide](./agents.md): driving gapit from agents
