# Phase-0 promoter robustness #2 — in-silico motif mutagenesis (causal)

Probe (routing_concat, promoter vs intergenic) trained on ORIGINAL train split,
scoring GC-preserving scrambles of val promoters. ΔP = drop in mean P(promoter)
vs the unperturbed window. `ctrl` = matched downstream scramble (the control).

val promoters scored: 2,784

## Overall

| variant | scrambled region | mean P(promoter) | ΔP vs original | % called promoter |
|---|---|---|---|---|
| original | (none) | 0.757 | +0.000 | 79.6% |
| m35 | −35 hexamer | 0.749 | +0.008 | 79.5% |
| m10 | −10 hexamer | 0.753 | +0.004 | 80.2% |
| m10m35 | −35 + −10 | 0.738 | +0.019 | 79.2% |
| core | core promoter (−40..−2) | 0.710 | +0.047 | 75.9% |
| ctrl | downstream +10..+48 (control) | 0.667 | +0.089 | 70.3% |

## Decision

- **Causal use of core promoter: NO** (ΔP core +0.047 vs ctrl +0.089; ≥ 2.0× control)
- Causal use of −35/−10 hexamers: NO (ΔP m10m35 +0.019 vs ctrl +0.089)

## Per-organism ΔP (specificity check)

| organism | n | original | m35 | m10 | m10m35 | core | ctrl |
|---|---|---|---|---|---|---|---|
| bsubtilis_168 | 133 | +0.000 | -0.019 | -0.002 | -0.021 | -0.030 | +0.046 |
| ecoli_K12_MG1655 | 1722 | +0.000 | +0.011 | +0.005 | +0.022 | +0.047 | +0.080 |
| hvolcanii_DS2 | 929 | +0.000 | +0.007 | +0.003 | +0.020 | +0.059 | +0.113 |

Reading: a large ΔP for `core`/`m10m35` with a near-zero ΔP for `ctrl` means
the model's promoter call causally depends on the promoter motif, not on a
co-occurring element elsewhere in the window.
