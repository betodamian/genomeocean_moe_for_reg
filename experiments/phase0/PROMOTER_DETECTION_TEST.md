# Confound-free promoter detection — design & rationale (2026-06-30)

After the robustness tests showed the promoter-vs-intergenic result was confounded by
downstream genic content, we designed a test that isolates promoter detection with **no
confound that can manufacture a positive result**. This documents the design analysis we
did *before* running anything.

## Why a literal "genomic-scan detection" test cannot be confound-free

A promoter is statistically entangled with its neighborhood (gene-proximal, AT-rich,
regulatory). Any **real** non-promoter negative differs from a promoter in ways correlated
with — but not identical to — the motif, so you can never match *everything*. The only way
to match everything is a **synthetic** negative (shuffle), which then introduces a
naturalness/out-of-distribution confound (the model detects "synthetic," not "no promoter").
That tension is fundamental, so we do not claim a confound-free *genomic-scan* metric.

## The test that IS confound-free (within-class dinucleotide-shuffle)

Reframe to the precise question that admits an airtight test: **does GO-MoE encode
promoter-specific structure that no composition, region, or naturalness confound explains?**

    signal = MCC(P vs P') − MCC(I vs I')

- **P** = real promoter window, **upstream-only** (−80..+1; `window_seq[70:151]`, TSS at idx 150) — no transcript.
- **P'** = **dinucleotide-preserving shuffle** of P (Altschul-Erikson; exact mono- + di-nucleotide composition, same length/endpoints; only the *positional* motif/structure destroyed).
- **I** = the same upstream slice of the intergenic decoys; **I'** = its dinuc-shuffle.

### Every confound is eliminated

| confound | why it's gone |
|---|---|
| downstream genic (the one we found) | upstream-only window — no transcript in any class |
| GC / mononucleotide | shuffle preserves exactly → invisible in P-vs-P' and I-vs-I' |
| dinucleotide structure | dinuc-preserving shuffle preserves exactly → invisible too |
| naturalness / OOD | a generic "shuffled-ness" cue appears equally in P-vs-P' and I-vs-I' → **cancels** in the subtraction |
| region-type entanglement | we **never compare P to I** — each class is compared to its own shuffle |
| train/val leakage | committed ≤60% cluster split; shuffles inherit their parent's split |

What survives the subtraction is **only** promoter-specific higher-order structure (the
−35/−10/spacer/UP architecture). And it answers the MoE question cleanly: if routing's
`(P−P') − (I−I')` exceeds embedding's, the promoter signal is **routing-carried**,
confound-free (P1).

### The single residual (minor, non-inflating)

If promoter regions had systematically different *sequence complexity* than intergenic for a
non-promoter reason, it could nudge the baseline — but complexity is largely mono/di
composition, which is preserved, so the effect is small and not systematically inflating. It
cannot manufacture a false positive.

### Scope (honest)

This tests "encodes promoter-specific structure beyond all controlled confounds" — the
strongest motif-isolating *detection* claim. It is **not** a genomic-scan deployment metric
(that needs real negatives and carries the unavoidable region entanglement). The two answer
different questions; this one is the airtight half.

### Generality

The identical `(X vs X') − (bg vs bg')` design applies to RBS (SD region) and Rho (rut
region) — a cleaner version of the earlier SD test — so it can only help the other elements.

## Implementation
- `pipeline/build_promoter_dinuc_test.py` — builds windows + verified dinuc-shuffles → `data/phase0/promoter_dinuc_ALL.tsv` (56,196 rows; 0 shuffle mismatches; GC identical real vs shuffle).
- `pipeline/phase0_extract_promoter_dinuc.py` + `submit_phase0_promoter_dinuc.sh` — frozen-MoE feature extraction (cluster).
- `pipeline/phase0_promoter_dinuc_analysis.py` — the MCC(P-P') − MCC(I-I') analysis, bootstrap CIs, routing-vs-embedding, per organism.

## RESULTS (2026-06-30, cluster job 23523691, frozen MoE)

MCC [95% bootstrap CI], committed ≤60% split:

| view | P vs P' | I vs I' (baseline) | promoter signal | clean CI? |
|---|---|---|---|---|
| embedding | 0.365 [0.340,0.388] | 0.319 [0.295,0.344] | +0.045 | no (overlap) |
| routing | 0.043 [0.018,0.069] | 0.057 [0.030,0.083] | −0.014 | no |
| **concat** | **0.388 [0.364,0.412]** | **0.336 [0.312,0.360]** | **+0.052** | **YES** |

Per-organism (concat signal): B. subtilis **+0.200**, H. volcanii +0.056, E. coli −0.027.

### Two clean conclusions

1. **There IS a genuine, confound-free promoter signal — but small.** concat +0.052 with
   non-overlapping CIs (0.364 > 0.360). So GO-MoE *does* encode promoter-specific structure
   beyond every controlled confound — it is not purely the downstream-genic artifact. Driven
   by B. subtilis (+0.200; canonical σA −35/−10, "bland" intergenic baseline 0.11); E. coli
   shows none because its intergenic sequence is itself highly structured (baseline 0.42),
   masking any promoter excess.

2. **The MoE routing channel does NOT carry it (P1 FAILS, confound-free).** routing signal
   −0.014 (P-vs-P' 0.043 ≈ its baseline 0.057 — both ~0). The promoter signal that exists
   lives in the **embedding / residual stream** (dense-accessible), not routing.

### This OVERTURNS the earlier confounded reading

The earlier zonal test (on the *confounded* full windows) tentatively read "routing 0.33 >
embedding 0.29 at the core → P1 holds." That was an artifact: routing is excellent at the
**genic-vs-intergenic structural axis** (Junho's specialty), which the confounded task
rewarded. Confound-free, routing is **blind to the promoter motif itself** — exactly the
pattern found for the short Shine-Dalgarno motif. So the honest, consistent story:

> The MoE router specializes for long structural classes (tRNA/rRNA/genic context) but does
> **not** allocate a routing channel to short regulatory motifs (SD, and now promoters under
> a fully confound-free test). Detection-wise, a small genuine promoter signal exists, but it
> is in the residual stream a dense model also has — not a MoE-necessity result for promoters.
