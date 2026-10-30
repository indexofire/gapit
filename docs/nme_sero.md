# Meningococcal serogroup typing (`nme_sero`)

`gapit` ships a fetch-on-demand **`nme_sero`** database for the
gene-presence part of *Neisseria meningitidis* typing:

- serogroup A / B / C / W / X / Y (W/Y when both amplicons are present)
- `ctrA` capsule marker
- `sodC` single-copy depth control

This matches gapit's positioning: **which known genes are present**. MLST/cc
is left to gMLST. Allele-level finetyping (PorA VR1/VR2, FetA, PorB,
BAST/MenDeVAR) is a separate sequence-typing problem and is deliberately
**not** part of this database; strain typing belongs to the gMLST database,
and antigen allele IDs to an external finetyping pipeline.

## Install

```bash
gapit db fetch nme_sero
gapit db list | grep nme_sero
gapit db info -d nme_sero
```

The database is downloaded at fetch time (GPL-3.0 / PubMLST terms); it is not
bundled inside the MIT wheel. The content snapshot used for validation is
upstream commit `35c4a5bff24fbc3625f92ce575d7c2e99ecd4378` (2026-02-19).

## Serogroup from contigs

```bash
gapit screen isolate.fna -d nme_sero
gapit typing isolate.tsv
```

`nme_sero/typing.json` declares the `men_serogroup` scheme. Its control gene
is `sodC`; an isolate without `sodC` is reported as untypeable rather than
silently assigned.

## Serogroup from FASTQ + sodC depth control

`nme_sero/screen.json` sets:

```json
"reads": {
  "min_identity": 80.0,
  "min_mapq": 20,
  "min_depth": 5.0,
  "min_relative_depth": 0.10,
  "relative_to": "sodC"
}
```

so a plain run gates presence on the gene:sodC depth ratio:

```bash
gapit screen -d nme_sero -1 R1.fastq.gz -2 R2.fastq.gz
```

`sodC` is a conserved single-copy Neisseria Cu-Zn superoxide dismutase used as
the relative depth denominator. It is not an absolute species-specific marker;
use `--min-relative-depth` / `--relative-to` to override the gate per run.

## Why no `gapit meningotype` subcommand?

`gapit typing` designates a phenotype from **gene presence/absence** rules.
PorA/FetA/PorB/BAST are **allele-level** calls: they need a protein allele
database, exact full-length peptide matching, and a BAST tuple lookup. That is
outside the generic `screen -> typing` contract, so it is not part of this
database. Use the gMLST database for strain typing and an external finetyping
pipeline for antigen allele IDs.

## Can `gapit typing` resolve the EX7E motif?

Not directly. `gapit typing` reads the screen table (gene, coverage, identity),
so it never sees the aligned codon and cannot translate the synG EX7E window.
The generic rule engine has no motif/allele primitive.

It can be approximated with **allele-specific probe records**: put the P, G and
S variants of a flanked probe into the database, raise the rule's
`identity_floor` to ~99.5, and let presence/absence choose the phenotype (this
is the "allele probe trick" used by the `typing` cookbook). That approach only
recognizes the exact probe sequences, is sensitive to silent variation in the
flank, and cannot report novel alleles; it is also fragile when a longer
serogroup amplicon overlaps the probe because abricate's `-culling_limit 1` can
shadow the probe hit.

The faithful implementation is upstream `menwy`: align `synG`, extract the
919–945 window (positions relative to the 3114 bp synG reference, **not** the
genome), translate it, and read residue 3 (P = W, G = Y, S = W/Y). That is
allele-level sequence interpretation and belongs in the dedicated finetyping
tool, not in this gene-presence database.

A W/Y call is therefore not silent about the gap: the rule's `notes` (JSON), the
`NOTES` column (TSV), and the Markdown table carry an explicit `NEXT STEP`
entry naming the EX7E window, the P/G/S mapping, and `menwy` as the follow-up
tool. The note also says the coordinates are synG-relative, because using
919–945 directly on a genome would be wrong.

An LLM can understand the instruction and route the sample to `menwy`, but it
should not try to "translate" the window without the aligned synG sequence or
the tool: the operation needs alignment, strand handling, and the bacterial
genetic code. Without `menwy` (or an equivalent allele caller), the correct
agent behavior is to stop and report that W vs Y is unresolved, not to guess
from the W/Y co-detection. The same guidance is available through the MCP
`typing` tool: pass the screen table text as `table` and read
`files[].phenotypes.men_serogroup.notes`.

## Caveats

- **W/Y** is called when W and Y amplicons are co-detected. The upstream synG
  EX7E motif refinement is not implemented; confirm rare W/Y variants with
  PCR or the upstream `menwy` tool. See the section above for why the generic
  `typing` engine is not the right place for it.
- `ctrA` is a capsule-locus marker, not a universal virulence gene; capsule-null
  isolates may lack it.
- `sodC` is a conserved single-copy control, not a formal species-specific
  marker.
- Database content is GPL-3.0 (meningotype) + PubMLST Neisseria terms; fetched
  on demand, not redistributed in the wheel.

## Validation

On the upstream A/B/C/W/X/Y test genomes, `gapit screen -d nme_sero` +
`gapit typing` returns A, B, C, W, X, Y respectively; `sodC` is present in all
six and `ctrA` in all six.
