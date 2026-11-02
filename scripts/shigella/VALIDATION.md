# Shigella serotype validation

## Datasets

- **S. boydii**: 144 assemblies (NCBI accessions from
  `~/data/shi/boydii/analysis/list`). References: ShigaTyper 2.0.5
  (`shigatyper_result.tsv`) and ShigaPass 1.5.0
  (`shigapass/ShigaPass_summary.csv`).
- **S. flexneri / S. sonnei / S. dysenteriae**: 669 assemblies from
  `~/data/shi/{flexneri,sonnei,dysenteriae}/01_fna` (312 + 293 + 64).
  The directory name is only a collection label, **not ground truth**: a few
  genomes collected under `flexneri/` are reference sonnei, and the accuracy
  tables are keyed by the independent reference call.
  Reference: ShigeiFinder 1.3.5 assembly calls
  (`~/data/shi/analysis/shigeifinder_all.tsv`, combining the 144
  existing boydii calls with the 669 fresh calls), with ShigaPass
  1.5.0 spot checks (`~/data/shi/analysis/shigapass_spot/`).
  ShigeiFinder leaves undetermined or reference-ambiguous genomes out of the
  accuracy denominator.

Commands:

```bash
gapit screen -F shi_all.fofn -d shi_sero -D /tmp/shi-dd -j 16 -o screen.tsv
gapit typing screen.tsv -D /tmp/shi-dd -o typing.tsv
python scripts/shigella/validate_serotyper.py \
  typing.tsv ~/data/shi/boydii/analysis/shigatyper_result.tsv
python scripts/shigella/validate_serotyper.py \
  typing.tsv ~/data/shi/boydii/analysis/shigapass/ShigaPass_summary.csv
python scripts/shigella/validate_serotyper.py \
  typing.tsv ~/data/shi/analysis/shigeifinder_all.tsv
```

The validator auto-detects the ShigaTyper two-column TSV, the semicolon
ShigaPass summary, and the ShigeiFinder `#SAMPLE`/`SEROTYPE` TSV.

## Collection directory is not truth

The four directories are collection folders, and their names do not reliably
identify the serogroup. ShigeiFinder's composition of the same 813 genomes:

| Collection directory | Genomes | ShigeiFinder reference species |
|---|---:|---|
| `boydii/` | 144 | boydii 90, flexneri 5, dysenteriae 6, undetermined 43 |
| `dysenteriae/` | 64 | dysenteriae 59, undetermined 5 |
| `flexneri/` | 312 | flexneri 291, sonnei 6, undetermined 15 |
| `sonnei/` | 293 | sonnei 291, undetermined 2 |

Examples: `GCA_001906365` sits in `boydii/` but is `SF2a`;
`GCA_001906425` is `SD3`; `GCF_042138705` sits in `flexneri/` but is `SS`.
Every accuracy table below is therefore keyed by the independent reference
call, never by the directory name.

## Results

All 813 genomes against the ShigeiFinder reference:

| Reference species | Determined | Exact | Ambiguous | Fallback | Mismatch |
|---|---:|---:|---:|---:|---:|
| S. boydii | 90 | 75 (83.3%) | 0 | 4 (4.4%) | 11 (12.2%) |
| S. flexneri | 296 | 293 (99.0%) | 0 | 0 | 3 (1.0%) |
| S. sonnei | 297 | 289 (97.3%) | 6 (2.0%) | 2 (0.7%) | 0 |
| S. dysenteriae | 65 | 64 (98.5%) | 0 | 1 (1.5%) | 0 |
| **All 813** | **748** | **721 (96.4%)** | 6 (0.8%) | 7 (0.9%) | 14 (1.9%) |

The boydii ShigeiFinder row is lower than the other species because the three
references disagree on a small set of boydii genomes, not because the marker
calls fail:

| S. boydii reference | Determined | Exact | Ambiguous | Fallback | Mismatch |
|---|---:|---:|---:|---:|---:|
| ShigaTyper 2.0.5 | 90 | 78 (86.7%) | 0 | 6 (6.7%) | 6 (6.7%) |
| ShigaPass 1.5.0 | 112 | 99 (88.4%) | 0 | 13 (11.6%) | 0 |
| ShigeiFinder 1.3.5 | 90 | 75 (83.3%) | 0 | 4 (4.4%) | 11 (12.2%) |

Of the 11 boydii/ShigeiFinder mismatches: 8 are SB6 vs SB10 (ShigaPass calls
all 8 SB10, and gapit calls SB10 from the intact `wbaM` 5' end), 2 are
SB1 vs SB20 (ShigaPass calls both SB20, and gapit calls SB20 from the
heparinase marker), and 1 is `SB13-atypical` vs SB13. The remaining
ShigeiFinder disagreements are not evidence that gapit's marker call is wrong.

## Remaining differences

### S. boydii SB6 vs SB10

ShigeiFinder calls eight genomes SB6; gapit calls SB10 from an intact 5'
`wbaM` alignment. ShigaPass also calls all eight SB10, matching the IS629
mechanism (Senchenkova et al. 2005). `GCA_003573425` is the 3'-partial case
(`wbaM 1-831/957`): its alignment starts at the 5' end, so it is SB10, not
an IS629 interruption.

### S. boydii SB1 vs SB20 (2)

`GCA_000815705` and `GCA_001905945` carry `Sb1_wzx` + `heparinase`; gapit
calls SB20, and ShigaPass also calls SB20. ShigeiFinder calls SB1 (it does
not use the heparinase distinction), so this is a reference-definition
difference rather than a gapit error.

### S. boydii 13-atypical

`GCA_001933075` is `SB13-atypical` in ShigeiFinder, `Not Shigella/EIEC` in
ShigaPass, and `S. boydii 13` in gapit (the control-exempt SB13 rule inherited
from ShigaTyper). This is an atypical/controversial serotype; the call is
reported as SB13 with the rule's note, not silently normalized to atypical.

### S. flexneri 2b -> 2a

`GCA_001906365.1` is labeled 2b by ShigaTyper, but the assembly contains no
`gtrX` at all (`blastn`, `blastx`, and `minimap2 asm20` against the 1266 bp
marker return nothing). ShigaTyper's own `report.json` evidence also lacks
`gtrX`; ShigeiFinder calls the genome `SF2a`. The 2a call stays; this is an
upstream label error.

### S. flexneri 2av vs 2a (3 genomes)

`GCF_056954245/485/565` carry `Xv` at ~99% identity with 70-86% coverage.
gapit calls the finer 2av; ShigeiFinder calls 2a, and ShigaPass reports
`SF1-5` with an `Unknown` POAC subtype. The 2av call is kept as a refinement.

### S. sonnei ambiguous (6)

These assemblies carry both S. flexneri (`Sf_wzx` + `gtrII`) and S. sonnei
(`Ss_methylase` + `Ss_wzx`) marker sets. ShigeiFinder calls them SS; gapit
keeps the explicit mixed-marker ambiguity rather than guessing.

### S. sonnei form II fallback (2)

`GCF_040493375` and `GCF_040495735` have `Ss_methylase` but no `cadA` call.
The form II rule keeps `cadA` because an S. dysenteriae 10 genome in the same
set also carries `Ss_methylase`; dropping the guard would turn that genome
ambiguous. ShigaPass spot checks call both SS; the conservative fallback
stays.

### Other fallbacks

The all-813 ShigeiFinder run has seven fallbacks: the two S. sonnei form II
cases above, four S. boydii (three weak-rfb SB4/SB8 assemblies that ShigaPass
still resolves, one ipaH-negative SB4), and one S. dysenteriae provisional
(the panel has no `SDprovBEDP02-5104` marker). These are marker-sensitivity
or control-limit cases, not mis-calls.

### S. dysenteriae 16 fallback (ShigaPass)

Two genomes (`GCF_002949345/555`) are called SD16 by ShigaPass from MLST +
fliC + CRISPR with only 2.8% rfb coverage. gapit has no Sd16 rule and its
`SdProv_wzx` marker does not match, so it reports unresolved. The other
ShigaPass fallbacks are similar weak-rfb cases; for example,
`GCF_000815515` is called SB1 from 12.2% rfb coverage, and gapit stays
unresolved.

### ipaH control

ShigaTyper's `ipaH_c` is one allele. Divergent Shigella ipaH alleles can lose
the full-gene coverage gate while retaining either the 5' or 3' part. The
bundle therefore splits `ipaH_c` into two non-overlapping presence fragments
(`1-640` and `641-780`). This recovered 16 S. flexneri fallbacks in the
non-boydii set without changing the 144-genome S. boydii result.

## Scope

This is an 813-genome smoke validation (144 S. boydii + 669 other Shigella),
not a replacement for the 4,879-genome ShigaPass validation or the
219-genome ShigaTyper study. ShigeiFinder is an independent assembly
reference, not a gold standard; every disagreement is listed above rather
than hidden in the aggregate.
