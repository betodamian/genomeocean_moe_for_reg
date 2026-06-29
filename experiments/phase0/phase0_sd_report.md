# Phase-0 RBS SD re-analysis — biophysical ΔG labels

Fixes the data defect in PHASE0_FINDINGS.md (circular regex labels + cross-organism GC leak). SAME frozen features; labels are now upstream:anti-SD hybridization ΔG (compute_sd_deltaG.py). WITHIN-organism evaluation so GC/taxonomy cannot leak. **P1 question: routing_only > embedding_only?**

## Macro-average across organisms

| view | (A) ΔG regression Spearman ρ | (B) strong-vs-weak SD MCC |
|---|---|---|
| gc_only | 0.017 | 0.021 |
| kmer4 | 0.279 | 0.266 |
| embedding_only | 0.136 | 0.117 |
| routing_only | 0.053 | 0.061 |
| routing_concat | 0.142 | 0.112 |

- **routing_only > embedding_only (P1):** regression NO (0.053 vs 0.136); binary NO (0.061 vs 0.117)
- **routing_concat beats no-model baselines (gc/kmer):** NO (ρ 0.142 vs gc 0.017, kmer 0.279)

## (A) ΔG regression — Spearman ρ per organism

| organism | n_val | gc_only | kmer4 | embedding_only | routing_only | routing_concat |
|---|---|---|---|---|---|---|
| bsubtilis_168 | 872 | 0.072 | 0.218 | 0.148 | 0.080 | 0.165 |
| ccrescentus_NA1000 | 770 | 0.004 | 0.326 | 0.226 | 0.015 | 0.217 |
| ecoli_BL21_DE3 | 610 | -0.030 | 0.369 | 0.054 | 0.038 | 0.071 |
| ecoli_K12_MG1655 | 421 | -0.002 | 0.353 | 0.101 | 0.052 | 0.093 |
| hvolcanii_DS2 | 307 | 0.026 | 0.189 | 0.080 | -0.014 | 0.067 |
| mtuberculosis_H37Rv | 710 | 0.049 | 0.355 | 0.245 | 0.071 | 0.246 |
| saureus_HG001 | 431 | -0.002 | 0.144 | 0.100 | 0.128 | 0.137 |

## (B) strong-vs-weak SD (per-organism ΔG terciles) — MCC per organism

| organism | n_val | ΔG cuts | gc_only | kmer4 | embedding_only | routing_only | routing_concat |
|---|---|---|---|---|---|---|---|
| bsubtilis_168 | 618 | [-11.5, -9.4] | 0.070 | 0.219 | 0.115 | 0.096 | 0.147 |
| ccrescentus_NA1000 | 503 | [-7.2, -4.5] | 0.005 | 0.340 | 0.216 | 0.027 | 0.172 |
| ecoli_BL21_DE3 | 419 | [-7.3, -5.2] | -0.011 | 0.306 | 0.074 | -0.010 | 0.047 |
| ecoli_K12_MG1655 | 280 | [-7.6, -5.6] | 0.073 | 0.371 | 0.014 | 0.146 | 0.079 |
| hvolcanii_DS2 | 213 | [-6.3, -3.9] | 0.010 | 0.117 | 0.033 | 0.033 | -0.002 |
| mtuberculosis_H37Rv | 475 | [-8.0, -4.2] | 0.020 | 0.323 | 0.234 | 0.026 | 0.256 |
| saureus_HG001 | 308 | [-10.2, -7.4] | -0.021 | 0.184 | 0.129 | 0.111 | 0.082 |
