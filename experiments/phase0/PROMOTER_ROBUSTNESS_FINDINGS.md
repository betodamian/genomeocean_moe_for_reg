# Phase-0 promoter robustness — findings (2026-06-30)

**Question (the "mystery element" check):** the Phase-0 promoter task (promoter vs
same-context intergenic decoy) scored routing_concat MCC 0.63 and passed the GC gate.
But is the model detecting the **promoter motif**, or a **co-occurring feature**? Two
frozen-model controls (one A40 job, `submit_phase0_promoter_robustness.sh`):

1. zonal per-position localization (`phase0_features_promoters_zonal.npz`, 0/28,098 token mismatches)
2. in-silico GC-preserving motif mutagenesis (`phase0_features_promoters_mutagenesis.npz`)

## Verdict: the promoter-vs-intergenic signal is substantially CONFOUNDED by downstream genic content

Both tests agree, and they identify the confound: a promoter window has a **gene/transcript
starting at the TSS** (its downstream half is genic), while the v1 decoy is **fully
intergenic**. So the model can separate them by *"is there a transcript downstream"* rather
than *"is there a −35/−10 promoter upstream."* This is the gene-vs-non-gene shortcut leaking
back through the **downstream half** of the window — the exact confound the same-context decoy
was meant to remove. (`build_promoter_windows.py` flagged the v1 decoy as provisional:
"the harder upstream-of-gene decoy is a later refinement.")

### Test #1 — zonal localization (val MCC, promoter vs decoy)

| zone | concat MCC | reading |
|---|---|---|
| full window | 0.629 | the headline number |
| **down (+6..+150, transcript body)** | **0.597** | **most of the signal is here** |
| farup (−150..−40) | 0.459 | |
| core promoter (−40..−1) | 0.432 | the actual promoter motif region |
| −10 / −35 / TSS | 0.34 / 0.32 / 0.34 | |

**Not localized to the promoter:** core 0.432 < down 0.597 (margin −0.164). The single
most discriminative zone is the **downstream transcript body**, not the −35/−10 motif.

### Test #2 — motif mutagenesis (ΔP = drop in mean P(promoter) when a region is scrambled)

| scrambled region | ΔP |
|---|---|
| −35 + −10 hexamers | +0.019 |
| core promoter (−40..−2) | +0.047 |
| **downstream control (+10..+48)** | **+0.089** |

**No causal dependence on the motif:** scrambling the promoter hexamers barely moves the
prediction (+0.019), while scrambling a *matched downstream control* moves it **more**
(+0.089). The trained probe leans on downstream/flanking content, not the promoter motif.
Same pattern in all three organisms (B. subtilis even shows ~0/negative ΔP for the motif
scrambles).

## What is and isn't true

- **Real but modest motif signal exists.** At the core-promoter zone alone (downstream
  excluded), concat MCC is 0.43 — far above the GC baseline (0.04) and the SD-region
  composition (kmer4 0.23) — and **routing 0.33 > embedding 0.29 (P1 holds at the motif).**
  So GenomeOcean-MoE *does* carry promoter-motif information in the routing channel; it is
  just much weaker than the 0.63 headline, and the full-window probe does not rely on it.
- **The headline 0.63 was inflated** by the downstream genic shortcut.
- **The robustness tests did exactly their job** — caught a confound *before* it became a
  false headline. This is the §4 / §12 "structural-shortcut as a confound to defeat" working.

## Scope — this likely affects the other "vs intergenic" tasks too

RBS-TIS and Rho are also feature-centered with genic content downstream of the site and
fully-intergenic decoys, so the **same downstream confound probably inflates
`rbs_TIS_vs_intergenic` and `rho_t1_vs_intergenic`**. The Tier-2 **same-position** task
`rho_t1_vs_intrinsic` (Rho vs intrinsic terminator, both in a 3′ context) is **immune** —
and notably it is the task where routing most dominated embedding (0.82 vs 0.41). That
contrast is the cleaner MoE-necessity headline.

## Recommended next step (decoy fix → re-run)

Rebuild the promoter (and RBS/Rho) **same-context decoys to match downstream genic content**
— e.g. windows that are upstream-of-a-gene / shifted into a transcript so both classes have
the same genic-vs-intergenic profile and differ only in the regulatory motif — then re-run
Phase-0. Until then, report the promoter result as **confounded**, with the de-confounded
core-zone estimate (concat ≈ 0.43, routing > embedding) as the honest motif-detection number.

Artifacts: `phase0_promoter_zonal_report.md`, `phase0_promoter_mutagenesis_report.md`
(+ JSONs). Features regenerable via `pipeline/submit_phase0_promoter_robustness.sh`.
