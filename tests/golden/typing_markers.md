---
schema: gapit.typing_result/1
tool: gapit 0.5.1
created_at: 2026-10-02T12:00:00Z
db: markers
source:
  - screen.tsv
files: 3
---

# gapit typing results

| FILE | SCHEME | PHENOTYPE | GENES | CONFIDENCE | SCORE | RUNNER_UP | NOTES |
|---|---|---|---|---|---|---|---|
| exact.fa | pathotype | EHEC | marker_a;marker_b;marker_c;marker_d | high | 1.0000 | EPEC (0.0000) |  |
| exact.fa | toxin | - | marker_a;marker_b;marker_c;marker_d | ambiguous | 1.0000 | - | ambiguous: Toxin1 (1.0000), Toxin2 (1.0000) |
| partial.fa | pathotype | EPEC | marker_a;marker_b | high | 1.0000 | EHEC (0.0000) |  |
| partial.fa | toxin | Toxin1 | marker_a;marker_b | high | 1.0000 | Toxin2 (0.0000) |  |
| mutated.fa | pathotype | EHEC | marker_a;marker_b;marker_c;marker_d | high | 1.0000 | EPEC (0.0000) |  |
| mutated.fa | toxin | Toxin2 | marker_a;marker_b;marker_c;marker_d | high | 1.0000 | Toxin1 (0.0000) |  |
