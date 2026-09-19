# Custom databases

`gapit db build` turns any FASTA of reference genes into a fully built gapit-native database in one command: `records.jsonl`, the `sequences` projection, the BLAST index, the manifest. This page is a worked tour with real session output. For providers, `db fetch`, and the datadir layout, see [Databases](./databases.md).

## When you need this

- You screen for genes no provider ships: a private marker panel, lab-curated variants.
- You want a small database with only the genes you report on, so outputs stay short.
- You have an abricate-era `~~~` FASTA and want it as a native gapit database.

## Command reference

```text
gapit db build NAME FASTA [OPTIONS]
```

| Option | Value | Meaning |
|---|---|---|
| `NAME` (argument) | text | Target database name (created under the datadir). |
| `FASTA` (argument) | path | Input FASTA: plain, abricate `~~~`, or `gapit\|` headers, detected per record (`.gz`/`.bz2` accepted). |
| `--tsv` | path | Metadata TSV: header row with `gene`/`accession`/`function` columns. |
| `--dbtype` | `nucl` or `prot` | Force the molecule type (default: auto-detect from the sequences). |
| `--datadir` | path | Database directory (default: `$GAPIT_DATADIR`, then `~/.local/share/gapit/db`). |
| `--description` | text | Default product for records whose FASTA header has no description text. |
| `--force` | flag | Overwrite the database if it already exists. |
| `--quiet` | flag | Silence stderr diagnostics. |

Build progress lines go to stderr, the one-line JSON receipt to stdout. The receipt carries `db`, `records`, `dbtype`, and `destination`, the same fields as `db fetch`. Screen blocks below show stdout only; `gapit screen` prints its `Processing:` lines to stderr (see [Screening](./screen.md)).

Every example below was executed in one session inside a scratch directory, so receipts quote its absolute path. To reproduce, run from the repo root:

```bash
export PATH="$PWD/.pixi/envs/default/bin:$PATH"
mkdir -p /tmp/opencode/customdb-docs && cd /tmp/opencode/customdb-docs
```

The input genes are synthetic, generated with Python's `random` module (a fixed seed keeps this page reproducible; never use real sequences you can't share):

```python
import random

rng = random.Random(11)
g1 = "".join(rng.choice("ACGT") for _ in range(240))
rng = random.Random(12)
g2 = "".join(rng.choice("ACGT") for _ in range(240))
open("genes.fa", "w").write(
    f">labcur1 synthetic tetracycline efflux pump\n{g1}\n"
    f">labcur2 synthetic macrolide esterase\n{g2}\n"
)
```

```console
$ grep '>' genes.fa
>labcur1 synthetic tetracycline efflux pump
>labcur2 synthetic macrolide esterase
```

## Worked examples

### 1. Plain FASTA to first hit

Build a database named `labgenes` into a local datadir, then screen a query contig that carries the first 200 bases of `labcur1` plus flanking bases:

```console
$ gapit db build labgenes genes.fa --datadir ./db
gapit: generated /tmp/opencode/customdb-docs/db/labgenes/sequences
gapit: self-check passed for labgenes
gapit: BLAST index built (nucl)
gapit: minimap2 index built
{"db":"labgenes","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labgenes"}
$ gapit screen query1.fa --db labgenes --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query1.fa	contigA	31	230	+	labcur1	1-200/240	=============..	0/0	83.33	100.00	labgenes		synthetic tetracycline efflux pump	
```

The hit row decodes the plain header: gene `labcur1`, product from the header description, empty `ACCESSION` and `RESISTANCE` because a plain header carries neither. `%COVERAGE` is 83.33 because the query holds 200 of the gene's 240 bases; the trailing `..` in the coverage map sketches those 40 uncovered bases.

### 2. Attaching metadata with --tsv

A metadata TSV adds accession and function classes per gene. Header row mandatory, `gene` column mandatory, `accession` and `function` optional per file:

```console
$ cat meta1.tsv
gene	accession	function
labcur1	LAB-0001	tetracycline
labcur2	LAB-0002	macrolide
$ gapit db build labmeta genes.fa --datadir ./db --tsv meta1.tsv
gapit: generated /tmp/opencode/customdb-docs/db/labmeta/sequences
gapit: self-check passed for labmeta
gapit: BLAST index built (nucl)
gapit: minimap2 index built
{"db":"labmeta","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labmeta"}
$ gapit screen query1.fa --db labmeta --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query1.fa	contigA	31	230	+	labcur1	1-200/240	=============..	0/0	83.33	100.00	labmeta	LAB-0001	synthetic tetracycline efflux pump	tetracycline
```

The `function` column is `;`-separated when a gene belongs to several classes. Both classes land in the `RESISTANCE` column as one string:

```console
$ cat meta2.tsv
gene	accession	function
labcur1	LAB-0001	virulence;marker
labcur2	LAB-0002	macrolide
$ gapit db build labmulti genes.fa --datadir ./db --tsv meta2.tsv
gapit: generated /tmp/opencode/customdb-docs/db/labmulti/sequences
gapit: self-check passed for labmulti
gapit: BLAST index built (nucl)
gapit: minimap2 index built
{"db":"labmulti","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labmulti"}
$ gapit screen query1.fa --db labmulti --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query1.fa	contigA	31	230	+	labcur1	1-200/240	=============..	0/0	83.33	100.00	labmulti	LAB-0001	synthetic tetracycline efflux pump	virulence;marker
```

### 3. abricate-style ~~~ headers

Files exported from abricate conventions build without changes. Fields split on `~~~`, and the description after the last field becomes the product:

```console
$ grep '>' legacy.fa
>oldlab~~~tetA_lab~~~SYN-100~~~TETRACYCLINE synthetic tetracycline pump
>oldlab~~~ermX~~~SYN-101~~~MACROLIDE synthetic methylase
$ gapit db build oldlab legacy.fa --datadir ./db
gapit: generated /tmp/opencode/customdb-docs/db/oldlab/sequences
gapit: self-check passed for oldlab
gapit: BLAST index built (nucl)
gapit: minimap2 index built
{"db":"oldlab","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/oldlab"}
$ gapit screen query_tet.fa --db oldlab --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query_tet.fa	contigB	1	240	+	tetA_lab	1-240/240	===============	0/0	100.00	100.00	oldlab	SYN-100	synthetic tetracycline pump	TETRACYCLINE
```

The first `~~~` field is the original database name and is discarded: `db` is always the NAME you are building. Gene, accession, and resistance all flow through to the screen output.

### 4. Roundtripping through gapit| headers

Native records are portable. Regenerate a FASTA from the `records.jsonl` written by example 1 (one `gapit|` header per record, sequence and product verbatim), then rebuild it as a new database:

```console
$ grep '>' ported.fa
>gapit|db=labgenes|gene=labcur1|acc=|func= synthetic tetracycline efflux pump
>gapit|db=labgenes|gene=labcur2|acc=|func= synthetic macrolide esterase
$ gapit db build labgenes_rt ported.fa --datadir ./db
gapit: generated /tmp/opencode/customdb-docs/db/labgenes_rt/sequences
gapit: self-check passed for labgenes_rt
gapit: BLAST index built (nucl)
gapit: minimap2 index built
{"db":"labgenes_rt","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labgenes_rt"}
$ gapit screen query1.fa --db labgenes_rt --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query1.fa	contigA	31	230	+	labcur1	1-200/240	=============..	0/0	83.33	100.00	labgenes_rt		synthetic tetracycline efflux pump	
```

The screen output is identical to example 1 (the `db` value in the header is retargeted to the new NAME). Empty `acc=` and `func=` segments are fine, and a `gapit|` header carrying accession and function values roundtrips them the same way.

### 5. Protein databases (blastx)

Amino acid input is detected automatically. Note the receipt says `"dbtype":"prot"` and the `minimap2 index built` step is gone (minimap2 indexes nucleotide only):

```console
$ grep '>' toxins.faa
>toxA synthetic pore-forming toxin
$ gapit db build toxins toxins.faa --datadir ./db
gapit: generated /tmp/opencode/customdb-docs/db/toxins/sequences
gapit: self-check passed for toxins
gapit: BLAST index built (prot)
{"db":"toxins","records":1,"dbtype":"prot","destination":"/tmp/opencode/customdb-docs/db/toxins"}
$ gapit db build toxins2 toxins.faa --datadir ./db --dbtype prot
gapit: generated /tmp/opencode/customdb-docs/db/toxins2/sequences
gapit: self-check passed for toxins2
gapit: BLAST index built (prot)
{"db":"toxins2","records":1,"dbtype":"prot","destination":"/tmp/opencode/customdb-docs/db/toxins2"}
```

Screening a protein database runs `blastx`, so the query file must be nucleotide. Here the query contig carries the toxin's coding sequence, back-translated one codon per amino acid, with flanking bases:

```console
$ gapit screen query_toxin.fa --db toxins --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query_toxin.fa	contigT	31	390	+	toxA	1-120/120	===============	0/0	100.00	100.00	toxins		synthetic pore-forming toxin
```

Reading the row: `START`/`END` are nucleotide coordinates on the contig (30 bp flank, then the 360 bp CDS), while `COVERAGE` and `%COVERAGE` count amino acids over the 120-residue protein. `--dbtype nucl|prot` exists for the edge case where the letter alphabet alone would mislead the heuristic; for ordinary inputs the auto-detect and the explicit flag agree.

### 6. Compressed input

`.gz` (and `.bz2`) inputs build exactly like plain ones:

```console
$ gapit db build labgz genes.fa.gz --datadir ./db
gapit: generated /tmp/opencode/customdb-docs/db/labgz/sequences
gapit: self-check passed for labgz
gapit: BLAST index built (nucl)
gapit: minimap2 index built
{"db":"labgz","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labgz"}
```

### 7. Iterating: rebuild with --force

Add a third gene to the FASTA and rebuild `labgenes` from example 1. An existing database is never overwritten silently:

```console
$ grep '>' genes_v2.fa
>labcur1 synthetic tetracycline efflux pump
>labcur2 synthetic macrolide esterase
>labcur3 synthetic fosfomycin thiolase
$ gapit db build labgenes genes_v2.fa --datadir ./db
{"schema":"gapit.error/1","code":"DB_ALREADY_EXISTS","message":"won't overwrite existing database labgenes (use --force)","context":{"db":"labgenes"}}
$ echo $?
4
```

The envelope is the standard `gapit.error/1` shape (exit 4, database error).

`--force` deletes and rebuilds in place. The metadata TSV from example 2 still applies:

```console
$ gapit db build labgenes genes_v2.fa --datadir ./db --tsv meta1.tsv --force
gapit: generated /tmp/opencode/customdb-docs/db/labgenes/sequences
gapit: self-check passed for labgenes
gapit: BLAST index built (nucl)
gapit: minimap2 index built
{"db":"labgenes","records":3,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labgenes"}
$ gapit screen query1.fa --db labgenes --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query1.fa	contigA	31	230	+	labcur1	1-200/240	=============..	0/0	83.33	100.00	labgenes	LAB-0001	synthetic tetracycline efflux pump	tetracycline
```

The receipt now reports `records: 3`. Edit `records.jsonl` or the FASTA, rebuild with `--force`, and every downstream artifact regenerates.

### 8. Metadata edge rules

Duplicate rows for one gene: the **first row wins**, a warning goes to stderr (`--quiet` silences it), and the build succeeds. The screen output proves which row survived:

```console
$ cat dup.tsv
gene	accession	function
labcur1	KEEP-1	tetracycline
labcur1	LOST-2	macrolide
$ gapit db build labdup genes.fa --datadir ./db --tsv dup.tsv
WARNING: duplicate gene 'labcur1' in metadata TSV: keeping the first row
gapit: generated /tmp/opencode/customdb-docs/db/labdup/sequences
gapit: self-check passed for labdup
gapit: BLAST index built (nucl)
gapit: minimap2 index built
{"db":"labdup","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labdup"}
$ gapit screen query1.fa --db labdup --datadir ./db
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
query1.fa	contigA	31	230	+	labcur1	1-200/240	=============..	0/0	83.33	100.00	labdup	KEEP-1	synthetic tetracycline efflux pump	tetracycline
```

A TSV row naming a gene the FASTA does not carry is warned about and skipped. The FASTA is the truth for what exists:

```console
$ cat ghost.tsv
gene	accession	function
labcur1	LAB-0001	tetracycline
ghostgene	G-999	virulence
$ gapit db build labghost genes.fa --datadir ./db --tsv ghost.tsv
WARNING: gene 'ghostgene' in metadata TSV not found in FASTA: skipped
gapit: generated /tmp/opencode/customdb-docs/db/labghost/sequences
gapit: self-check passed for labghost
gapit: BLAST index built (nucl)
gapit: minimap2 index built
{"db":"labghost","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labghost"}
$ grep -c '"gene"' db/labghost/records.jsonl
2
```

### 9. Where builds land without --datadir

The datadir resolves per call: `--datadir`, then `$GAPIT_DATADIR`, then `~/.local/share/gapit/db`. This example points `GAPIT_DATADIR` at a scratch directory (unset, the same command writes into the default `~/.local/share/gapit/db`):

```console
$ export GAPIT_DATADIR=/tmp/opencode/customdb-default/db
$ gapit db build labenv genes.fa
gapit: generated /tmp/opencode/customdb-default/db/labenv/sequences
gapit: self-check passed for labenv
gapit: BLAST index built (nucl)
gapit: minimap2 index built
{"db":"labenv","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-default/db/labenv"}
$ gapit list
DATABASE	SEQUENCES	DBTYPE	DATE
labenv	2	nucl	2026-Sep-19
```

Unlike screening, a build creates a missing datadir instead of failing, so a fresh machine bootstraps on the first build. The scratch datadir was deleted after this capture.

### 10. Finishing the pipeline: two samples, one matrix

Screen two contig files against `labmeta` from example 2, save the report tables, and fold them into a summary matrix with [gapit summary](./summary.md):

```console
$ gapit screen sample1.fa --db labmeta --datadir ./db > sample1.tsv
$ gapit screen sample2.fa --db labmeta --datadir ./db > sample2.tsv
$ cat sample1.tsv
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
sample1.fa	SAM001	41	280	+	labcur1	1-240/240	===============	0/0	100.00	100.00	labmeta	LAB-0001	synthetic tetracycline efflux pump	tetracycline
$ cat sample2.tsv
#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
sample2.fa	SAM002	16	255	+	labcur2	1-240/240	===============	0/0	100.00	100.00	labmeta	LAB-0002	synthetic macrolide esterase	macrolide
$ gapit summary sample1.tsv sample2.tsv
#FILE	NUM_FOUND	labcur1	labcur2
sample1.tsv	1	100.00	.
sample2.tsv	1	.	100.00
```

## Rules and edge behavior

### Header auto-detection

The header kind is detected per record, so one FASTA can mix all three formats:

| FASTA header | Gene | Accession | Function | Product |
|---|---|---|---|---|
| `>mixplain plain header gene` | `mixplain` | none | none | `plain header gene` |
| `>oldlab~~~fromtilde~~~SYN-100~~~TETRACYCLINE tilde gene` | `fromtilde` | `SYN-100` | `TETRACYCLINE` | `tilde gene` |
| `>gapit\|db=elsewhere\|gene=fromtag\|acc=SYN-200\|func=ampicillin;gentamicin tagged gene` | `fromtag` | `SYN-200` | `ampicillin`, `gentamicin` | `tagged gene` |
| `>secondplain` (bare) | `secondplain` | none | none | `--description`, else `secondplain` |

For a plain header, the text after the id is the product when present. When absent, `--description TEXT` fills in, and with neither the gene name becomes the product:

```console
$ grep '>' bare.fa
>baregene
$ gapit db build bare_desc bare.fa --datadir ./db --description "lab-curated reference gene"
gapit: generated /tmp/opencode/customdb-docs/db/bare_desc/sequences
gapit: self-check passed for bare_desc
gapit: BLAST index built (nucl)
gapit: minimap2 index built
{"db":"bare_desc","records":1,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/bare_desc"}
$ grep -o '"product":"[^"]*"' db/bare_desc/records.jsonl
"product":"lab-curated reference gene"
$ gapit db build bare_plain bare.fa --datadir ./db --force
gapit: generated /tmp/opencode/customdb-docs/db/bare_plain/sequences
gapit: self-check passed for bare_plain
gapit: BLAST index built (nucl)
gapit: minimap2 index built
{"db":"bare_plain","records":1,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/bare_plain"}
$ grep -o '"product":"[^"]*"' db/bare_plain/records.jsonl
"product":"baregene"
```

A header that starts with `gapit|` but is not valid tagged syntax fails the build with exit 4:

```console
$ gapit db build labbroken broken_tag.fa --datadir ./db
{"schema":"gapit.error/1","code":"HEADER_MALFORMED","message":"malformed gapit/v1 sequence header: segment_without_key","context":{"seqid":"gapit|nonsense","reason":"segment_without_key"}}
```

The mixed file from the table above builds each record by its own rules, and the generated `sequences` projection shows the retargeting (`db` is always the NAME built, the original `elsewhere` tag is dropped):

```console
$ gapit db build labmix mixed.fa --datadir ./db 2>/dev/null
{"db":"labmix","records":4,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labmix"}
$ grep '>' db/labmix/sequences
>gapit|db=labmix|gene=mixplain|acc=|func= plain header gene
>gapit|db=labmix|gene=fromtilde|acc=SYN-100|func=TETRACYCLINE tilde gene
>gapit|db=labmix|gene=fromtag|acc=SYN-200|func=ampicillin%3Bgentamicin tagged gene
>gapit|db=labmix|gene=secondplain|acc=|func= secondplain
```

Records keep input order, so `records.jsonl` line N is the Nth FASTA record. Duplicate gene names are allowed; a FASTA carrying `twingene` twice builds a two-record database:

```console
$ gapit db build labtwin dupname.fa --datadir ./db --force
gapit: generated /tmp/opencode/customdb-docs/db/labtwin/sequences
gapit: self-check passed for labtwin
gapit: BLAST index built (nucl)
gapit: minimap2 index built
{"db":"labtwin","records":2,"dbtype":"nucl","destination":"/tmp/opencode/customdb-docs/db/labtwin"}
$ grep -c '"gene":"twingene"' db/labtwin/records.jsonl
2
```

### dbtype resolution

`--dbtype` wins when given; otherwise the abricate mol-type heuristic decides from the input sequences themselves (`nucl` unless the letters say protein). Protein databases get a `.pin` BLAST index and no minimap2 index, and screening goes through `blastx` with a nucleotide query (example 5).

### What lands in the database directory

```console
$ ls db/labmeta
gapit-manifest.json
records.jsonl
sequences
sequences.mmi
sequences.ndb
sequences.nhr
sequences.nin
sequences.njs
sequences.not
sequences.nsq
sequences.ntf
sequences.nto
```

`records.jsonl` is the editable truth and `sequences` its `gapit/v1`-header projection; the `sequences.n*` files are the BLAST index and `sequences.mmi` the minimap2 index (nucleotide only). Field-level detail for both file contracts is in [Databases](./databases.md). Sequences are stored verbatim, no provider-style normalization: this is your curated truth, and a manifest with `source_urls: ["local"]` certifies the build.

## Troubleshooting

| Symptom | Error (exit code) | Fix |
|---|---|---|
| Build refuses an existing database | `DB_ALREADY_EXISTS` (4) | Rebuild with `--force` (example 7). |
| Metadata TSV rejected | `METADATA_MALFORMED` (5) | The first row must be a header containing a `gene` column. |
| Build fails on an empty input | `BUILD_INVALID` (4) | The FASTA parsed to zero records; check the file has `>` headers. |
| FASTA or TSV path missing | `INPUT_NOT_FOUND` (5) | The envelope's `context.file` names the unreadable path. |
| Build fails on a `gapit\|` header | `HEADER_MALFORMED` (4) | Fix the tagged header; plain and `~~~` headers can't use `gapit\|` syntax. |
| Screen against the custom db finds nothing | none, exit 0 | Empty result tables are normal output. Loosen `--minid`/`--mincov` or check the query matches the strand; see [Screening](./screen.md). |

Error envelopes are documented in [Outputs](./outputs.md); the full screening flag set in [Screening](./screen.md).
