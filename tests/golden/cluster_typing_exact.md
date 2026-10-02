---
schema: gapit.cluster/1
tool: gapit 0.5.0
created_at: 2026-09-30T12:00:00Z
db: cps
preset: asm20
min_gene_cov: 90.0
min_gene_id: 90.0
min_cluster_cov: 96.0
threads: 1
---

# gapit cluster screening report

## `exact.fa`

| File | Best locus | Type | Phenotype | Coverage% | Identity% | Present | Partial | Missing |
|---|---|---|---|---|---|---|---|---|
| exact.fa | locusA | KL101 | K101 | 100.00 | 100.00 | 3 | 0 | - |

Phenotype `K101` (score 1.0000, high confidence)

Best locus `locusA` (K antigen locus A, type KL101):

| Gene | Start | End | Strand | Coverage% | Identity% | Verdict |
|---|---|---|---|---|---|---|
| wzx | 101 | 500 | - | 100.00 | 100.00 | present |
| gtrA | 501 | 900 | + | 100.00 | 100.00 | present |
| manC | 901 | 1200 | - | 100.00 | 100.00 | present |

