# Phase-0 promoter robustness — ruling out a co-occurring "mystery element"

**Question.** The Phase-0 smoke test showed the frozen MoE separates promoter windows
from same-context intergenic decoys (routing_concat MCC 0.63, GC-matched gate PASS). But
"separates promoter windows" is not yet "detects the promoter": a feature that merely
*co-occurs* with promoters — some other element prevalent in promoter windows — could be
driving the classifier. These two cheap tests (frozen model, no retraining) localize the
signal in space and test its causal dependence on the promoter motif.

The same-context decoy and GC-matched baseline already exclude the trivial confounds
(gene-vs-non-gene, GC). These tests target an *unknown co-located* confound.

---

## Test #1 — Zonal (per-position) localization  `phase0_extract_promoter_zonal.py`

Capture **per-token** routing + embedding and pool them within zones of the 300-bp
window (TSS at index 150; token→bp coords from BPE strings, the method that gave
0/20,696 mismatches on RBS). Train a promoter-vs-decoy probe per zone, score val MCC.

| zone | bp index | meaning |
|---|---|---|
| core | 110–150 | core promoter (−40..−1): bacterial −35/−10 **and** archaeal TATA/BRE |
| minus35 / minus10 | 112–124 / 136–146 | bacterial hexamers |
| tss | 146–156 | transcription start |
| down | 156–300 | transcript body — **control** (no promoter motif) |
| farup | 0–110 | far upstream — **control** |
| full | all | reproduces the original whole-window pool |

**Decision:** if `MCC(core) ≫ MCC(down), MCC(farup)` the discriminative signal sits **at
the promoter region**, so a mystery element elsewhere in the window is not what drives the
call. Also reports `routing_core > emb_core` (P1: the MoE channel carries the motif).

## Test #2 — In-silico motif mutagenesis (causal)  `phase0_extract_promoter_mutagenesis.py`

The strongest test. GC-**preservingly** scramble the promoter motif in each val promoter
(exact base composition kept → any drop is not a GC artifact), re-run the frozen forward
pass, and score the **already-trained** promoter probe. Variants: `m35`, `m10`, `m10m35`,
`core` (−40..−2), and `ctrl` (matched downstream +10..+48 — the control).

**Decision:** if `ΔP(core)` and `ΔP(m10m35)` are ≥ 2× `ΔP(ctrl)`, the promoter call
**causally depends on the promoter motif**, not on a co-occurring element. Per-organism
specificity check: bacteria should respond to `m10`/`m35`; archaeal *H. volcanii* to
`core` (its TATA/BRE architecture differs) — if the bacterial hexamer scrambles hit only
bacteria, that itself confirms motif-specific causality.

---

## How to run

1. Sync `pipeline/` to the cluster and submit (one A40 job, ~10 min):
   ```
   sbatch pipeline/submit_phase0_promoter_robustness.sh
   ```
2. Pull the two NPZs back to `data/datasets/phase0/` (tar | ssh; scp is broken), then:
   ```
   .venv/bin/python pipeline/phase0_promoter_zonal_analysis.py \
       --zonal data/datasets/phase0/phase0_features_promoters_zonal.npz \
       --all data/datasets/promoters/ALL.tsv --out experiments/phase0
   .venv/bin/python pipeline/phase0_promoter_mutagenesis_analysis.py \
       --orig data/datasets/phase0/phase0_features_promoters.npz \
       --mut  data/datasets/phase0/phase0_features_promoters_mutagenesis.npz \
       --out experiments/phase0
   ```
   Reports: `phase0_promoter_zonal_report.md`, `phase0_promoter_mutagenesis_report.md`.

## What each outcome means

- **Both pass** (signal localized to core + scrambling the motif collapses the call):
  strong evidence the model detects the promoter itself — a mystery-element confound is
  ruled out. This hardens the headline promoter result before the Week-3 baselines.
- **Localized but mutagenesis-insensitive:** the signal is at the promoter region but the
  model leans on context around the motif rather than the hexamers — report honestly.
- **Not localized / control scramble drops as much:** a co-occurring feature *is*
  contributing — flag it and design a matched-confound control before claiming detection.

These extend, and feed into, the formal §9b structural partial-out and the promoter-vs-RBS
acid test planned for Weeks 3–4.
