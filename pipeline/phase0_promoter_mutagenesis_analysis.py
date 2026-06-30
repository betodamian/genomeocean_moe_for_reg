#!/usr/bin/env python3
"""
Analyze promoter MUTAGENESIS features (robustness test #2) — the causal control.

Train the Phase-0 promoter probe (routing_concat, promoter vs intergenic decoy) on the
ORIGINAL features' train split, then score each GC-preserving motif-scramble variant of
the val promoters. If destroying the promoter motif collapses P(promoter) far more than
scrambling a matched downstream control, the model is CAUSALLY using the promoter — not
a co-occurring "mystery" element.

ΔP(variant) = mean P(promoter | original) − mean P(promoter | variant), over val positives.

Decision:
  ΔP(core) and ΔP(m10m35) ≥ FACTOR × ΔP(ctrl)  → causal use of the promoter motif.
  Per-organism: bacteria should respond to m10/m35; archaeal H. volcanii to core (its
  TATA/BRE architecture differs) — a built-in specificity check.

Output: <out>/phase0_promoter_mutagenesis_report.md + ..._results.json
Run: .venv/bin/python pipeline/phase0_promoter_mutagenesis_analysis.py \
        --orig data/datasets/phase0/phase0_features_promoters.npz \
        --mut  data/datasets/phase0/phase0_features_promoters_mutagenesis.npz \
        --out experiments/phase0
"""
import argparse, json, os
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

VARIANTS = ["original", "m35", "m10", "m10m35", "core", "ctrl"]
FACTOR = 2.0   # motif scramble must drop P at least 2x more than the control scramble


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orig", required=True, help="phase0_features_promoters.npz (probe training)")
    ap.add_argument("--mut", required=True, help="phase0_features_promoters_mutagenesis.npz")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    O = np.load(args.orig, allow_pickle=False)
    label = O["label"].astype(str); split = O["split"].astype(str)
    y = (label == "positive").astype(int)
    tr = split == "train"
    Xtr = np.concatenate([O["emb"], O["routing"]], axis=1)[tr]
    sc = StandardScaler(); Xtr_s = sc.fit_transform(Xtr)
    clf = LogisticRegression(max_iter=5000, class_weight="balanced", n_jobs=-1, C=1.0)
    clf.fit(Xtr_s, y[tr])

    M = np.load(args.mut, allow_pickle=False)
    org = M["organism"].astype(str)
    def pscore(v):
        X = np.concatenate([M[f"emb_{v}"], M[f"routing_{v}"]], axis=1)
        return clf.predict_proba(sc.transform(X))[:, 1]
    P = {v: pscore(v) for v in VARIANTS}

    base = P["original"]
    overall = {}
    for v in VARIANTS:
        dP = float(np.mean(base - P[v]))
        called = float(np.mean(P[v] >= 0.5))
        overall[v] = {"mean_P_promoter": round(float(np.mean(P[v])), 4),
                      "delta_P_vs_original": round(dP, 4),
                      "frac_called_promoter": round(called, 4)}

    dP_ctrl = overall["ctrl"]["delta_P_vs_original"]
    causal_core = overall["core"]["delta_P_vs_original"] >= FACTOR * max(dP_ctrl, 1e-3)
    causal_hex = overall["m10m35"]["delta_P_vs_original"] >= FACTOR * max(dP_ctrl, 1e-3)

    per_org = {}
    for o in sorted(set(org)):
        om = org == o
        if om.sum() < 20: continue
        b = base[om]
        per_org[o] = {v: round(float(np.mean(b - P[v][om])), 4) for v in VARIANTS}
        per_org[o]["n"] = int(om.sum())

    results = {"overall": overall, "per_organism": per_org,
               "decision": {"causal_core_promoter": bool(causal_core),
                            "causal_hexamers_m10m35": bool(causal_hex),
                            "delta_P_core": overall["core"]["delta_P_vs_original"],
                            "delta_P_m10m35": overall["m10m35"]["delta_P_vs_original"],
                            "delta_P_ctrl": dP_ctrl, "factor": FACTOR}}
    with open(os.path.join(args.out, "phase0_promoter_mutagenesis_results.json"), "w") as fh:
        json.dump(results, fh, indent=2)

    L = ["# Phase-0 promoter robustness #2 — in-silico motif mutagenesis (causal)\n",
         "\nProbe (routing_concat, promoter vs intergenic) trained on ORIGINAL train split,\n"
         "scoring GC-preserving scrambles of val promoters. ΔP = drop in mean P(promoter)\n"
         "vs the unperturbed window. `ctrl` = matched downstream scramble (the control).\n",
         f"\nval promoters scored: {len(base):,}\n",
         "\n## Overall\n\n| variant | scrambled region | mean P(promoter) | ΔP vs original | % called promoter |\n",
         "|---|---|---|---|---|\n"]
    region_desc = {"original": "(none)", "m35": "−35 hexamer", "m10": "−10 hexamer",
                   "m10m35": "−35 + −10", "core": "core promoter (−40..−2)",
                   "ctrl": "downstream +10..+48 (control)"}
    for v in VARIANTS:
        o = overall[v]
        L.append(f"| {v} | {region_desc[v]} | {o['mean_P_promoter']:.3f} | "
                 f"{o['delta_P_vs_original']:+.3f} | {o['frac_called_promoter']*100:.1f}% |\n")
    L.append("\n## Decision\n\n")
    L.append(f"- **Causal use of core promoter: {'YES' if causal_core else 'NO'}** "
             f"(ΔP core {overall['core']['delta_P_vs_original']:+.3f} vs ctrl {dP_ctrl:+.3f}; "
             f"≥ {FACTOR}× control)\n")
    L.append(f"- Causal use of −35/−10 hexamers: {'YES' if causal_hex else 'NO'} "
             f"(ΔP m10m35 {overall['m10m35']['delta_P_vs_original']:+.3f} vs ctrl {dP_ctrl:+.3f})\n")
    L.append("\n## Per-organism ΔP (specificity check)\n\n")
    L.append("| organism | n | " + " | ".join(VARIANTS) + " |\n")
    L.append("|---|---|" + "|".join(["---"]*len(VARIANTS)) + "|\n")
    for o, r in per_org.items():
        L.append(f"| {o} | {r['n']} | " + " | ".join(f"{r[v]:+.3f}" for v in VARIANTS) + " |\n")
    L.append("\nReading: a large ΔP for `core`/`m10m35` with a near-zero ΔP for `ctrl` means\n"
             "the model's promoter call causally depends on the promoter motif, not on a\n"
             "co-occurring element elsewhere in the window.\n")
    with open(os.path.join(args.out, "phase0_promoter_mutagenesis_report.md"), "w") as fh:
        fh.write("".join(L))
    print("".join(L))
    print(f"wrote {args.out}/phase0_promoter_mutagenesis_report.md + .json")


if __name__ == "__main__":
    main()
