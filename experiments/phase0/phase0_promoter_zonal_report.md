# Phase-0 promoter robustness #1 — zonal (per-position) localization

Logistic probe (promoter vs intergenic decoy) trained per ZONE on the committed
train split; val MCC. Localizes WHERE in the 300-bp window the discriminative
signal lives — to rule out a co-occurring confound outside the promoter.

val n=5594 (2784 promoters)

## val MCC by zone × feature view

| zone (bp) | gc | kmer4 | embedding | routing | concat |
|---|---|---|---|---|---|
| full (all) | -0.024 | 0.408 | 0.440 | 0.537 | 0.629 |
| core (110–150) | 0.042 | 0.227 | 0.292 | 0.332 | 0.432 |
| minus35 (112–124) | 0.086 | 0.251 | 0.248 | 0.215 | 0.321 |
| minus10 (136–146) | 0.083 | 0.246 | 0.265 | 0.219 | 0.344 |
| tss (146–156) | 0.043 | 0.141 | 0.276 | 0.221 | 0.343 |
| down (156–300) | 0.050 | 0.407 | 0.485 | 0.448 | 0.597 |
| farup (0–110) | 0.012 | 0.281 | 0.238 | 0.424 | 0.459 |

## Decision

- **Localized to core promoter: NO** (concat MCC core 0.432 − max(down 0.597, farup 0.459) = -0.164; margin ≥ 0.1)
- P1 at core (routing > embedding): YES (0.332 vs 0.292)

Reading: if `core` ≫ `down`/`farup`, the promoter call is driven by the
−35/−10/TSS region, not by a co-occurring element elsewhere in the window.
