#!/usr/bin/env python3
"""
Analyze the confound-free promoter-detection test (within-class dinucleotide-shuffle).

    signal = MCC(P vs P') - MCC(I vs I')        per feature view

P-vs-P' = can the model tell a real promoter from its dinucleotide-shuffle (GC + dinuc
identical; only the positional motif/structure destroyed)? I-vs-I' = the same for
intergenic = the generic "ordered-structure / naturalness" baseline. The DIFFERENCE is
promoter-specific structure that no composition / region / naturalness confound explains.

MoE-necessity (P1), confound-free: routing's signal > embedding's signal?

Bootstrap 95% CI on each MCC; the promoter signal is "real" if the P-vs-P' MCC clears
the I-vs-I' MCC with non-overlapping CIs.

Output: <out>/phase0_promoter_dinuc_report.md + ..._results.json
Run: .venv/bin/python pipeline/phase0_promoter_dinuc_analysis.py \
        --npz data/datasets/phase0/phase0_features_promoter_dinuc.npz --out experiments/phase0
"""
import argparse, json, os
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import matthews_corrcoef

RNG = np.random.default_rng(20260630)
N_BOOT = 1000
VIEWS = ["embedding", "routing", "concat"]


def build(view, emb, rout):
    if view == "embedding": return emb
    if view == "routing":   return rout
    return np.concatenate([emb, rout], axis=1)


def boot_mcc(yte, ypred):
    n = len(yte); vals = []
    for _ in range(N_BOOT):
        idx = RNG.integers(0, n, n)
        if len(np.unique(yte[idx])) < 2: continue
        vals.append(matthews_corrcoef(yte[idx], ypred[idx]))
    return (round(float(np.percentile(vals, 2.5)), 4), round(float(np.percentile(vals, 97.5)), 4))


def contrast(real_klass, shuf_klass, Z, view):
    """MCC distinguishing real vs its dinuc-shuffle, on the committed split."""
    klass = Z["klass"]; split = Z["split"]
    m = np.isin(klass, [real_klass, shuf_klass])
    X = build(view, Z["emb"], Z["routing"])[m]
    y = (klass[m] == real_klass).astype(int)
    sp = split[m]
    tr, te = sp == "train", sp == "val"
    sc = StandardScaler(); Xtr = sc.fit_transform(X[tr]); Xte = sc.transform(X[te])
    clf = LogisticRegression(max_iter=5000, class_weight="balanced", C=1.0, n_jobs=-1)
    clf.fit(Xtr, y[tr]); yp = clf.predict(Xte)
    mcc = round(float(matthews_corrcoef(y[te], yp)), 4)
    lo, hi = boot_mcc(y[te], yp)
    return {"mcc": mcc, "ci": [lo, hi], "n_val": int(te.sum())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    d = np.load(args.npz, allow_pickle=False)
    Z = {k: d[k] for k in d.files}
    Z["klass"] = Z["klass"].astype(str); Z["split"] = Z["split"].astype(str)
    Z["organism"] = Z["organism"].astype(str)

    res = {"views": {}}
    for v in VIEWS:
        P = contrast("P", "P_shuf", Z, v)        # promoter vs its shuffle
        I = contrast("I", "I_shuf", Z, v)        # intergenic baseline
        signal = round(P["mcc"] - I["mcc"], 4)
        clean = P["ci"][0] > I["ci"][1]          # non-overlapping CIs
        res["views"][v] = {"P_vs_Pshuf": P, "I_vs_Ishuf": I,
                           "promoter_specific_signal": signal, "clean_CI_separation": bool(clean)}

    # MoE-necessity, confound-free: routing signal vs embedding signal
    sig = {v: res["views"][v]["promoter_specific_signal"] for v in VIEWS}
    res["moe_necessity"] = {"routing_signal": sig["routing"], "embedding_signal": sig["embedding"],
                            "routing_gt_embedding": sig["routing"] > sig["embedding"]}

    # per-organism (concat view)
    per_org = {}
    for o in sorted(set(Z["organism"])):
        sub = {k: (Z[k][Z["organism"] == o] if k in ("klass", "split") else Z[k][Z["organism"] == o])
               for k in ("klass", "split")}
        Zo = {"emb": Z["emb"][Z["organism"] == o], "routing": Z["routing"][Z["organism"] == o],
              "klass": sub["klass"], "split": sub["split"]}
        try:
            P = contrast("P", "P_shuf", Zo, "concat"); I = contrast("I", "I_shuf", Zo, "concat")
            per_org[o] = {"P_vs_Pshuf": P["mcc"], "I_vs_Ishuf": I["mcc"],
                          "signal": round(P["mcc"] - I["mcc"], 4), "n_val": P["n_val"]}
        except Exception:
            pass
    res["per_organism_concat"] = per_org

    with open(os.path.join(args.out, "phase0_promoter_dinuc_results.json"), "w") as fh:
        json.dump(res, fh, indent=2)

    L = ["# Phase-0 confound-free promoter detection — within-class dinucleotide-shuffle\n",
         "\n`signal = MCC(P vs P') − MCC(I vs I')`, where X' is the dinucleotide-preserving\n"
         "shuffle of X (GC + dinucleotide composition identical; only positional motif/structure\n"
         "destroyed), on upstream-only windows (−80..+1). The intergenic baseline cancels the\n"
         "generic ordered-structure/naturalness cue; the difference is promoter-specific structure\n"
         "that no composition / region / naturalness confound can explain.\n",
         "\n## Result (MCC [95% CI])\n\n",
         "| view | P vs P' | I vs I' (baseline) | promoter signal | clean? |\n",
         "|---|---|---|---|---|\n"]
    for v in VIEWS:
        r = res["views"][v]; P, I = r["P_vs_Pshuf"], r["I_vs_Ishuf"]
        L.append(f"| {v} | {P['mcc']:.3f} [{P['ci'][0]:.3f},{P['ci'][1]:.3f}] | "
                 f"{I['mcc']:.3f} [{I['ci'][0]:.3f},{I['ci'][1]:.3f}] | "
                 f"**{r['promoter_specific_signal']:+.3f}** | {'YES' if r['clean_CI_separation'] else 'no'} |\n")
    mo = res["moe_necessity"]
    L.append(f"\n## MoE-necessity (confound-free P1)\n\n")
    L.append(f"- routing promoter-signal {mo['routing_signal']:+.3f} vs embedding {mo['embedding_signal']:+.3f} "
             f"→ routing > embedding: **{'YES' if mo['routing_gt_embedding'] else 'NO'}**\n")
    L.append("\n## Per-organism (concat)\n\n| organism | P vs P' | I vs I' | signal | n_val |\n|---|---|---|---|---|\n")
    for o, r in per_org.items():
        L.append(f"| {o} | {r['P_vs_Pshuf']:.3f} | {r['I_vs_Ishuf']:.3f} | {r['signal']:+.3f} | {r['n_val']} |\n")
    L.append("\nReading: a positive `promoter signal` with clean CI separation = the model encodes\n"
             "promoter-specific structure beyond all controlled confounds. routing>embedding = that\n"
             "structure is carried by the MoE routing channel (MoE-necessity), confound-free.\n")
    with open(os.path.join(args.out, "phase0_promoter_dinuc_report.md"), "w") as fh:
        fh.write("".join(L))
    print("".join(L))
    print(f"wrote {args.out}/phase0_promoter_dinuc_report.md + .json")


if __name__ == "__main__":
    main()
