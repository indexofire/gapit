---
schema: gapit.cluster/2
tool: gapit 0.5.6
created_at: 2026-10-01T12:00:00Z
db: vp_ok
preset: asm20
min_gene_cov: 90.0
min_gene_id: 90.0
min_cluster_cov: 96.0
threads: 1
---

# gapit cluster screening report

## `vp_o3k6.fa`

| FILE | SEROTYPE | O_GROUP | O_LOCUS | O_COVERAGE | O_IDENTITY | O_PRESENT | O_PARTIAL | O_MISSING | K_GROUP | K_LOCUS | K_COVERAGE | K_IDENTITY | K_PRESENT | K_PARTIAL | K_MISSING |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| vp_o3k6.fa | O3/O13:K6 | O3/O13 | OL3_or_O13 | 100.00 | 100.00 | 3 | 0 | - | K6 | KL6 | 100.00 | 100.00 | 3 | 0 | - |

Phenotype `o_group` = O3/O13 (score 1.0000, high confidence; winner locus `OL3_or_O13` coverage 100.00%, identity 100.00%, 3 present, 0 partial, missing -)
Phenotype `k_group` = K6 (score 1.0000, high confidence; winner locus `KL6` coverage 100.00%, identity 100.00%, 3 present, 0 partial, missing -)
Phenotype `serotype` = O3/O13:K6 (score 1.0000, high confidence)

Best locus `OL3_or_O13` (OL3, type O3/O13):

| Gene | Start | End | Strand | Coverage% | Identity% | Verdict |
|---|---|---|---|---|---|---|
| o3_wzx | 101 | 533 | + | 100.00 | 100.00 | present |
| o3_wzy | 601 | 1033 | - | 100.00 | 100.00 | present |
| o3_gmd | 1101 | 1400 | + | 100.00 | 100.00 | present |

