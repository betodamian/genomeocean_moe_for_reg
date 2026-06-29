# Phase-0 smoke test — ceiling benchmark (research_plan §4, P0)
Frozen GenomeOcean-MoE (pilot4 lr_match_dropout05). Linear probe on the committed 80/20 cluster-level split. **Gate:** `routing_concat` MCC lower-CI > `gc_only` MCC upper-CI (clears the GC-matched baseline with non-overlapping 95% CIs).

## Verdict: 4/5 tasks PASS the GC-matched gate

| task | n_val (pos) | routing_concat MCC | gc_only MCC | gate |
|---|---|---|---|---|
| promoter_vs_intergenic | 5594 (2784) | 0.629 [0.610,0.649] | -0.024 [-0.051,0.002] | **PASS** |
| rbs_TIS_vs_intergenic | 8261 (4121) | 0.807 [0.793,0.819] | 0.017 [-0.005,0.038] | **PASS** |
| rbs_SD_vs_UNSD | 4121 (2592) | 0.335 [0.305,0.364] | 0.280 [0.252,0.310] | FAIL |
| rho_t1_vs_intergenic | 1307 (179) | 0.698 [0.644,0.751] | 0.549 [0.499,0.601] | **PASS** |
| rho_t1_vs_intrinsic_ecoli | 318 (30) | 0.779 [0.658,0.890] | 0.096 [-0.014,0.196] | **PASS** |

## promoter_vs_intergenic

- element: `promoters`  positives=['positive']  negatives=['decoy']  restrict_org=None
- train: 22504 (11265 pos)  val: 5594 (2784 pos)
- gate: PASS  (routing_concat − gc_only MCC = +0.653)

| feature view | MCC [95% CI] | AUPRC | AUROC | F1 | acc |
|---|---|---|---|---|---|
| gc_only | -0.024 [-0.051,0.002] | 0.555 | 0.529 | 0.486 | 0.488 |
| kmer4 | 0.408 [0.383,0.432] | 0.749 | 0.775 | 0.700 | 0.704 |
| embedding_only | 0.440 [0.418,0.465] | 0.790 | 0.793 | 0.713 | 0.720 |
| routing_only | 0.537 [0.515,0.557] | 0.877 | 0.855 | 0.760 | 0.768 |
| routing_concat | 0.629 [0.610,0.649] | 0.910 | 0.900 | 0.809 | 0.814 |

routing_concat val MCC by organism: bsubtilis_168=0.513(n=271), ecoli_K12_MG1655=0.634(n=3445), hvolcanii_DS2=0.638(n=1878)

## rbs_TIS_vs_intergenic

- element: `rbs`  positives=['positive_SD', 'positive_UNSD']  negatives=['decoy_intergenic']  restrict_org=None
- train: 33131 (16575 pos)  val: 8261 (4121 pos)
- gate: PASS  (routing_concat − gc_only MCC = +0.790)

| feature view | MCC [95% CI] | AUPRC | AUROC | F1 | acc |
|---|---|---|---|---|---|
| gc_only | 0.017 [-0.005,0.038] | 0.523 | 0.537 | 0.501 | 0.508 |
| kmer4 | 0.463 [0.443,0.482] | 0.786 | 0.806 | 0.732 | 0.731 |
| embedding_only | 0.614 [0.597,0.630] | 0.859 | 0.879 | 0.808 | 0.807 |
| routing_only | 0.710 [0.695,0.726] | 0.943 | 0.938 | 0.855 | 0.855 |
| routing_concat | 0.807 [0.793,0.819] | 0.970 | 0.968 | 0.903 | 0.903 |

routing_concat val MCC by organism: bsubtilis_168=0.820(n=1741), ccrescentus_NA1000=0.878(n=1540), ecoli_BL21_DE3=0.768(n=1240), ecoli_K12_MG1655=0.799(n=854), hvolcanii_DS2=0.682(n=615), mtuberculosis_H37Rv=0.751(n=1412), saureus_HG001=0.898(n=859)

## rbs_SD_vs_UNSD

- element: `rbs`  positives=['positive_SD']  negatives=['positive_UNSD']  restrict_org=None
- train: 16575 (10274 pos)  val: 4121 (2592 pos)
- gate: FAIL  (routing_concat − gc_only MCC = +0.055)

| feature view | MCC [95% CI] | AUPRC | AUROC | F1 | acc |
|---|---|---|---|---|---|
| gc_only | 0.280 [0.252,0.310] | 0.738 | 0.666 | 0.697 | 0.646 |
| kmer4 | 0.394 [0.365,0.422] | 0.838 | 0.768 | 0.751 | 0.704 |
| embedding_only | 0.337 [0.311,0.368] | 0.804 | 0.729 | 0.725 | 0.675 |
| routing_only | 0.270 [0.243,0.298] | 0.761 | 0.669 | 0.683 | 0.636 |
| routing_concat | 0.335 [0.305,0.364] | 0.806 | 0.728 | 0.723 | 0.674 |

routing_concat val MCC by organism: bsubtilis_168=0.176(n=872), ccrescentus_NA1000=0.159(n=770), ecoli_BL21_DE3=0.193(n=610), ecoli_K12_MG1655=0.228(n=421), hvolcanii_DS2=0.130(n=307), mtuberculosis_H37Rv=0.263(n=710), saureus_HG001=0.148(n=431)

## rho_t1_vs_intergenic

- element: `rho`  positives=['positive_rho_t1', 'positive_rho_t1_rsr']  negatives=['decoy_intergenic']  restrict_org=None
- train: 5268 (714 pos)  val: 1307 (179 pos)
- gate: PASS  (routing_concat − gc_only MCC = +0.149)

| feature view | MCC [95% CI] | AUPRC | AUROC | F1 | acc |
|---|---|---|---|---|---|
| gc_only | 0.549 [0.499,0.601] | 0.568 | 0.901 | 0.595 | 0.836 |
| kmer4 | 0.540 [0.491,0.589] | 0.583 | 0.895 | 0.589 | 0.836 |
| embedding_only | 0.482 [0.414,0.536] | 0.539 | 0.845 | 0.553 | 0.839 |
| routing_only | 0.672 [0.624,0.720] | 0.845 | 0.969 | 0.707 | 0.897 |
| routing_concat | 0.698 [0.644,0.751] | 0.816 | 0.953 | 0.740 | 0.922 |

routing_concat val MCC by organism: ecoli_K12_MG1655=0.475(n=901), mtuberculosis_H37Rv=0.552(n=297)

## rho_t1_vs_intrinsic_ecoli

- element: `rho`  positives=['positive_rho_t1', 'positive_rho_t1_rsr']  negatives=['decoy_intrinsic']  restrict_org=ecoli_K12_MG1655
- train: 1224 (121 pos)  val: 318 (30 pos)
- gate: PASS  (routing_concat − gc_only MCC = +0.683)

| feature view | MCC [95% CI] | AUPRC | AUROC | F1 | acc |
|---|---|---|---|---|---|
| gc_only | 0.096 [-0.014,0.196] | 0.314 | 0.663 | 0.206 | 0.541 |
| kmer4 | 0.146 [0.008,0.277] | 0.239 | 0.710 | 0.245 | 0.767 |
| embedding_only | 0.411 [0.235,0.562] | 0.374 | 0.824 | 0.467 | 0.899 |
| routing_only | 0.815 [0.695,0.906] | 0.928 | 0.991 | 0.831 | 0.965 |
| routing_concat | 0.779 [0.658,0.890] | 0.878 | 0.969 | 0.800 | 0.962 |

routing_concat val MCC by organism: ecoli_K12_MG1655=0.779(n=318)
