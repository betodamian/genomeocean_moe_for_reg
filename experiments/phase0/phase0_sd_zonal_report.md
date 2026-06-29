# Phase-0 RBS SD — zonal (per-position) re-analysis

Within-organism ΔG regression (Spearman ρ). Routing/embedding pooled over the
SD-zone tokens ONLY (no whole-window averaging), vs the original full-window pool.

## Macro-average across organisms (Spearman ρ predicting SD ΔG)

| view | ρ |
|---|---|
| gc_sd | 0.209 |
| kmer4_sd | 0.870 |
| emb_full | 0.149 |
| routing_full | 0.058 |
| emb_sd | 0.397 |
| routing_sd | 0.013 |
| routing_start | 0.003 |
| concat_sd | 0.380 |

**Decision: (b) non-specialization — un-pooling does NOT recover SD; router does not encode it**

- routing_sd 0.013 vs routing_full 0.058 → dilution NO
- P1 at SD zone (routing_sd > emb_sd): NO (0.013 vs 0.397)
- concat_sd > kmer4_sd (model beats SD-region composition): NO (0.380 vs 0.870)

## Per-organism ρ

| organism | n_val | gc_sd | kmer4_sd | emb_full | routing_full | emb_sd | routing_sd | routing_start | concat_sd |
|---|---|---|---|---|---|---|---|---|---|
| bsubtilis_168 | 872 | 0.272 | 0.791 | 0.125 | 0.112 | 0.400 | 0.097 | 0.066 | 0.399 |
| ccrescentus_NA1000 | 770 | 0.198 | 0.910 | 0.246 | 0.019 | 0.516 | -0.069 | 0.032 | 0.491 |
| ecoli_BL21_DE3 | 610 | 0.234 | 0.902 | 0.085 | 0.051 | 0.346 | 0.043 | -0.018 | 0.327 |
| ecoli_K12_MG1655 | 421 | 0.287 | 0.875 | 0.105 | 0.061 | 0.274 | 0.080 | 0.058 | 0.264 |
| hvolcanii_DS2 | 307 | 0.102 | 0.813 | 0.105 | -0.025 | 0.311 | -0.064 | -0.056 | 0.283 |
| mtuberculosis_H37Rv | 710 | 0.089 | 0.910 | 0.245 | 0.111 | 0.535 | 0.022 | 0.015 | 0.530 |
| saureus_HG001 | 431 | 0.278 | 0.892 | 0.134 | 0.077 | 0.392 | -0.015 | -0.078 | 0.369 |
