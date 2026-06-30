# Phase-0 confound-free promoter detection — within-class dinucleotide-shuffle

`signal = MCC(P vs P') − MCC(I vs I')`, where X' is the dinucleotide-preserving
shuffle of X (GC + dinucleotide composition identical; only positional motif/structure
destroyed), on upstream-only windows (−80..+1). The intergenic baseline cancels the
generic ordered-structure/naturalness cue; the difference is promoter-specific structure
that no composition / region / naturalness confound can explain.

## Result (MCC [95% CI])

| view | P vs P' | I vs I' (baseline) | promoter signal | clean? |
|---|---|---|---|---|
| embedding | 0.365 [0.340,0.388] | 0.319 [0.295,0.344] | **+0.045** | no |
| routing | 0.043 [0.018,0.069] | 0.057 [0.030,0.083] | **-0.014** | no |
| concat | 0.388 [0.364,0.412] | 0.336 [0.312,0.360] | **+0.052** | YES |

## MoE-necessity (confound-free P1)

- routing promoter-signal -0.014 vs embedding +0.045 → routing > embedding: **NO**

## Per-organism (concat)

| organism | P vs P' | I vs I' | signal | n_val |
|---|---|---|---|---|
| bsubtilis_168 | 0.309 | 0.109 | +0.200 | 266 |
| ecoli_K12_MG1655 | 0.394 | 0.421 | -0.027 | 3444 |
| hvolcanii_DS2 | 0.431 | 0.375 | +0.056 | 1858 |

Reading: a positive `promoter signal` with clean CI separation = the model encodes
promoter-specific structure beyond all controlled confounds. routing>embedding = that
structure is carried by the MoE routing channel (MoE-necessity), confound-free.
