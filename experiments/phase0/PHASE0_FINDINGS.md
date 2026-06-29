# Phase-0 ceiling benchmark — findings & go/no-go (research_plan §4)

**Date:** 2026-06-29
**Model:** frozen GenomeOcean-MoE, pilot4 `lr_match_dropout05` stage-1 step-49999
(8 experts, top-k 2, 12 layers × 768-d; the same checkpoint behind Junho's routing
results). No fine-tuning — frozen forward pass + linear probe only.
**Features:** 84,310 windows (300 bp) → `embedding_only` (768-d last hidden, mean-pooled)
and `routing` (96-d mean router softmax over 12×8). `routing_concat` = concat (864-d).
Baselines: `gc_only` (1-d GC fraction), `kmer4` (256-d 4-mer freqs).
**Probe:** logistic regression on the committed 80/20 ≤60%-identity cluster split.
**Gate (P0):** `routing_concat` clears the GC-matched baseline with non-overlapping
95% bootstrap CIs on MCC.

---

## Verdict: **GO** — 4/5 same-context tasks pass; the 1 failure is the pre-registered RBS risk

| task | routing_concat MCC | gc_only MCC | routing_only | embedding_only | gate |
|---|---|---|---|---|---|
| promoter vs intergenic | **0.629** [0.610,0.649] | −0.024 | 0.537 | 0.440 | **PASS** |
| RBS-TIS vs intergenic | **0.807** [0.793,0.819] | 0.017 | 0.710 | 0.614 | **PASS** |
| RBS SD vs UNSD | 0.335 [0.305,0.364] | 0.280 | 0.270 | 0.337 | **FAIL** |
| Rho vs intergenic | **0.698** [0.644,0.751] | 0.549 | 0.672 | 0.482 | **PASS** |
| Rho vs intrinsic (E. coli) | **0.779** [0.658,0.890] | 0.096 | 0.815 | 0.411 | **PASS** |

The frozen model already contains strong, GC-independent discriminative signal for
promoters, translation-initiation sites, and Rho terminators — including the hardest
same-context contrast (Rho-dependent vs intrinsic terminators, both terminators in the
same 3′ context). Building the full annotator (§7–§9) is justified.

---

## P1 evidence (routing > embedding), seen here already at Phase 0

The plan's central MoE-necessity prediction (P1: `routing_concat` > `embedding_only`,
because a dense model has no routing channel) shows up in **4/5 tasks** before any of
the dedicated §9 analysis:

- promoter:        routing_only 0.537 **>** embedding_only 0.440 ; concat lifts to 0.629
- RBS-TIS:         routing_only 0.710 **>** embedding_only 0.614 ; concat 0.807
- Rho-intergenic:  routing_only 0.672 **>** embedding_only 0.482 ; concat 0.698
- Rho-intrinsic:   routing_only **0.815** **>** embedding_only 0.411 ; concat 0.779

The Rho-vs-intrinsic case is the most striking: distinguishing two kinds of terminator
in the same context is carried almost entirely by the **routing channel** (0.815),
while the residual-stream embedding a dense model would expose is far weaker (0.411).
This is the qualitative shape P1/P3 predict; it will be tested formally in Week 3–4
(paired folds, expert ablation), but the direction is already present.

---

## The one failure is exactly the pre-registered RBS risk (not a kill)

`RBS SD vs UNSD` — discriminating Shine-Dalgarno-led starts from unleadered (UNSD)
starts — fails the gate: routing_concat 0.335 barely exceeds gc_only 0.280 and the CIs
overlap. Two diagnostic facts:

1. **`kmer4` (0.394) beats every model view** here — the weak separability that exists
   is plain sequence composition, and the pooled model representation does not add to it.
2. **`gc_only` is already 0.280** — SD vs UNSD carries real GC/composition confound.

This is precisely the failure the plan flagged in advance (§5d): *"RBS is the highest
Phase-0 risk element — the SD core is 6–8 bp ≈ 1–2 tokens inside a 300 bp ≈ 60-token
window, so the pooled `routing_concat` can dilute the SD signal below detectability."*
Crucially, **RBS itself is not dead**: the model localizes translation-initiation sites
overwhelmingly (RBS-TIS MCC 0.807, per-organism 0.68–0.90 across all 7 genomes). It is
specifically the 6–8 bp SD-motif sub-signal that the 300-bp **pooled** view washes out.

**Pre-registered response (§4 failure-mode 2 → change windowing/feature views):** run
the RBS window-size + pooling sweep — 300/120/60 bp windows, and **per-position /
per-expert-bin** routing views instead of the mean-pooled fingerprint — to recover the
SD motif. This is a targeted refinement of one sub-task, not a Phase-0 failure.

---

## Per-organism robustness (confound check)

`routing_concat` validation MCC is consistent across organisms, so the wins are not an
artifact of one genome dominating:

- **promoter:** E. coli 0.634, H. volcanii 0.638, B. subtilis 0.513 (incl. the archaeal
  domain — promoter signal transfers across the bacteria↔archaea boundary).
- **RBS-TIS:** all 7 organisms 0.68–0.90 (S. aureus 0.898, C. crescentus 0.878).
- **Rho-intergenic:** E. coli 0.475, MTB 0.552 (both in-vivo phyla).

---

## Addendum (2026-06-29) — RBS SD label fix: data defect confirmed, then a deeper finding

The `SD vs UNSD` failure was first diagnosed as a **data defect** and fixed:
- *Circular labels:* the original SD/UNSD call was a regex (`sd_label()` = "upstream
  contains a GGAGG-like substring"), so the task was substring rediscovery — a k-mer
  counter beat the model by construction.
- *Cross-organism GC leak:* pooled, "SD" coincided with "low-GC species" (B. subtilis
  85% SD @ 40% GC vs H. volcanii 24% SD @ 63% GC); within-organism the GC gap was ±0.005.

**Fix applied** (`compute_sd_deltaG.py`, `phase0_sd_reanalysis.py`): relabel SD strength
as the biophysical hybridization ΔG between each gene's upstream and that organism's own
16S 3′ anti-SD (ViennaRNA), evaluated **within organism**. The continuous ΔG is
biologically correct (B. subtilis strongest median −10.3; archaeal H. volcanii / MTB
weakest −4.8 / −6.2). Same frozen features, new labels — no GPU re-run.

**Result (macro-avg over 7 organisms, Spearman ρ predicting ΔG):**

| view | ΔG regression ρ | strong-vs-weak SD MCC |
|---|---|---|
| gc_only | 0.017 | 0.021 |
| kmer4 | **0.279** | **0.266** |
| embedding_only | 0.136 | 0.117 |
| routing_only | 0.053 | 0.061 |
| routing_concat | 0.142 | 0.112 |

Two things are now clear and honest:
1. **The fix worked:** `gc_only` collapses to ~0 — the confound and circularity are gone.
2. **The deeper finding:** even with clean labels, the frozen **pooled** MoE does not
   expose SD — `kmer4` beats every model view, and **routing_only < embedding_only (P1
   fails for SD)**. The SD signal is real and lives in local composition (k-mers recover
   it), but mean-pooling over ~60 tokens washes out the 1–2-token SD motif, and the
   router carries even less of it than the residual stream a dense model already has.

This splits the remaining question in two, to be settled by the pre-registered method
fix: **(a) pooling dilution** (real SD signal, lost to mean-pooling → recoverable with
per-position / per-expert-bin routing + 60-bp windows) vs **(b) genuine non-
specialization** (the router simply does not allocate a channel to a short regulatory
motif, unlike Junho's structural tRNA/rRNA experts — a real negative for P3 on RBS-SD).
Either outcome is publishable; the per-position re-extraction distinguishes them.

## Next steps

1. **RBS SD method fix (pre-registered, decisive):** re-extract RBS at 60/120 bp windows
   + **per-position / per-expert-bin** routing views (not mean-pooled) and re-run the
   ΔG re-analysis. Distinguishes pooling-dilution (a) from non-specialization (b) above.
2. **Week 3 — detection vs baselines (P2):** run BPROM / Prodigal-RBS / RhoTermPredict
   and the dense gLMs (NT, DNABERT-2, Evo 2, ProkBERT, GO-dense) on the *same* held-out
   windows; compare AUPRC + boundary-F1.
3. **Week 4 — MoE necessity (P3/P4) & generalization (P8):** expert-ablation DiD on the
   passing elements; leave-one-organism-out cross-genome folds (already committed in
   `folds_loo.tsv`) for the unseen-split test.

**Artifacts:** `phase0_report.md` (full metric tables), `phase0_results.json` (machine-
readable). Features cached as `phase0_features_<element>.npz` (gitignored; regenerate via
`pipeline/submit_phase0.sh`).
