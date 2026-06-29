#!/usr/bin/env python3
"""
Phase-0 RBS SD re-analysis with biophysical ΔG labels (fixes the data defect found
in PHASE0_FINDINGS.md: the original SD-vs-UNSD labels were a circular regex and were
cross-organism GC-confounded).

Uses the SAME frozen features already extracted (phase0_features_rbs.npz) — only the
LABELS change — joined by window_id to the ΔG from compute_sd_deltaG.py. No GPU re-run.

Two WITHIN-ORGANISM evaluations (so the cross-organism GC/taxonomy leak cannot help),
per feature view {gc_only, kmer4, embedding_only, routing_only, routing_concat}:

  (A) REGRESSION (primary, threshold-free): predict continuous SD ΔG; report Spearman ρ
      on val. Directly asks how much SD-pairing information each representation carries.
  (B) BINARY strong-vs-weak SD: per-organism ΔG terciles (top third = strong SD,
      bottom third = weak/none, middle dropped) -> a clean, balanced, GC-controlled
      contrast; report MCC on val.

Headline question (P1): does `routing_only` carry SD information beyond `embedding_only`
(what a dense model also has), and do either beat the gc_only / kmer4 no-model baselines?

Output: <out>/phase0_sd_report.md + phase0_sd_results.json
Run: .venv/bin/python pipeline/phase0_sd_reanalysis.py \
        --features data/datasets/phase0/phase0_features_rbs.npz \
        --deltaG data/datasets/rbs/sd_deltaG.tsv --out experiments/phase0
"""
import argparse, csv, json, os
import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import matthews_corrcoef

VIEWS = ["gc_only", "kmer4", "embedding_only", "routing_only", "routing_concat"]
MIN_N = 40   # min train/val examples to evaluate an organism


def build_view(d, view):
    if view == "gc_only":        return d["gc"].reshape(-1, 1)
    if view == "kmer4":          return d["kmer4"]
    if view == "embedding_only": return d["emb"]
    if view == "routing_only":   return d["routing"]
    if view == "routing_concat": return np.concatenate([d["emb"], d["routing"]], axis=1)
    raise ValueError(view)


def load_deltaG(path):
    dg = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            dg[r["window_id"]] = float(r["dG"])
    return dg


def regression(Xtr, ytr, Xte, yte):
    sc = StandardScaler(); Xtr = sc.fit_transform(Xtr); Xte = sc.transform(Xte)
    m = Ridge(alpha=1.0).fit(Xtr, ytr)
    rho = spearmanr(yte, m.predict(Xte)).correlation
    return round(float(rho), 4)


def binary(Xtr, ytr, Xte, yte):
    sc = StandardScaler(); Xtr = sc.fit_transform(Xtr); Xte = sc.transform(Xte)
    clf = LogisticRegression(max_iter=5000, class_weight="balanced", n_jobs=-1).fit(Xtr, ytr)
    return round(float(matthews_corrcoef(yte, clf.predict(Xte))), 4)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", required=True)
    ap.add_argument("--deltaG", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    d = np.load(args.features, allow_pickle=False)
    data = {k: d[k] for k in d.files}
    ids = data["ids"].astype(str); org = data["organism"].astype(str)
    split = data["split"].astype(str)
    dgmap = load_deltaG(args.deltaG)

    # align ΔG to feature rows (positives only have ΔG)
    dG = np.array([dgmap.get(i, np.nan) for i in ids], dtype=np.float64)
    has = ~np.isnan(dG)

    organisms = sorted(set(org[has]))
    results = {"regression": {}, "binary": {}, "organisms": organisms}

    for o in organisms:
        om = has & (org == o)
        tr = om & (split == "train"); te = om & (split == "val")
        if tr.sum() < MIN_N or te.sum() < MIN_N:
            continue

        # (A) regression on continuous ΔG
        reg = {}
        for v in VIEWS:
            X = build_view(data, v)
            reg[v] = regression(X[tr], dG[tr], X[te], dG[te])
        results["regression"][o] = {"n_train": int(tr.sum()), "n_val": int(te.sum()), **reg}

        # (B) binary strong-vs-weak by per-organism terciles (fit on train ΔG)
        lo, hi = np.quantile(dG[tr], [1/3, 2/3])
        def lab(mask):
            y = np.full(mask.sum(), -1)
            dv = dG[mask]
            y[dv <= lo] = 0      # weak / no SD
            y[dv >= hi] = 1      # strong SD
            return y
        ytr_all, yte_all = lab(tr), lab(te)
        ktr, kte = ytr_all >= 0, yte_all >= 0
        if ktr.sum() < MIN_N or kte.sum() < MIN_N:
            continue
        tr_idx = np.where(tr)[0][ktr]; te_idx = np.where(te)[0][kte]
        ytr, yte = ytr_all[ktr], yte_all[kte]
        bina = {}
        for v in VIEWS:
            X = build_view(data, v)
            bina[v] = binary(X[tr_idx], ytr, X[te_idx], yte)
        results["binary"][o] = {
            "n_train": int(ktr.sum()), "n_val": int(kte.sum()),
            "dG_tercile_cuts": [round(float(lo), 2), round(float(hi), 2)], **bina}

    # macro-averages across organisms
    def macro(section):
        out = {}
        for v in VIEWS:
            vals = [results[section][o][v] for o in results[section]
                    if not np.isnan(results[section][o][v])]
            out[v] = round(float(np.mean(vals)), 4) if vals else float("nan")
        return out
    results["regression_macro"] = macro("regression")
    results["binary_macro"] = macro("binary")

    with open(os.path.join(args.out, "phase0_sd_results.json"), "w") as fh:
        json.dump(results, fh, indent=2)

    # ── report ──────────────────────────────────────────────────────────────
    L = ["# Phase-0 RBS SD re-analysis — biophysical ΔG labels\n",
         "\nFixes the data defect in PHASE0_FINDINGS.md (circular regex labels + "
         "cross-organism GC leak). SAME frozen features; labels are now upstream:anti-SD "
         "hybridization ΔG (compute_sd_deltaG.py). WITHIN-organism evaluation so GC/"
         "taxonomy cannot leak. **P1 question: routing_only > embedding_only?**\n"]

    rm, bm = results["regression_macro"], results["binary_macro"]
    L.append("\n## Macro-average across organisms\n\n")
    L.append("| view | (A) ΔG regression Spearman ρ | (B) strong-vs-weak SD MCC |\n")
    L.append("|---|---|---|\n")
    for v in VIEWS:
        L.append(f"| {v} | {rm[v]:.3f} | {bm[v]:.3f} |\n")
    p1_reg = rm["routing_only"] > rm["embedding_only"]
    p1_bin = bm["routing_only"] > bm["embedding_only"]
    beats_base = (rm["routing_concat"] > max(rm["gc_only"], rm["kmer4"]))
    L.append(f"\n- **routing_only > embedding_only (P1):** regression "
             f"{'YES' if p1_reg else 'NO'} ({rm['routing_only']:.3f} vs {rm['embedding_only']:.3f}); "
             f"binary {'YES' if p1_bin else 'NO'} ({bm['routing_only']:.3f} vs {bm['embedding_only']:.3f})\n")
    L.append(f"- **routing_concat beats no-model baselines (gc/kmer):** "
             f"{'YES' if beats_base else 'NO'} "
             f"(ρ {rm['routing_concat']:.3f} vs gc {rm['gc_only']:.3f}, kmer {rm['kmer4']:.3f})\n")

    L.append("\n## (A) ΔG regression — Spearman ρ per organism\n\n")
    L.append("| organism | n_val | " + " | ".join(VIEWS) + " |\n")
    L.append("|---|---|" + "|".join(["---"]*len(VIEWS)) + "|\n")
    for o in results["regression"]:
        r = results["regression"][o]
        L.append(f"| {o} | {r['n_val']} | " +
                 " | ".join(f"{r[v]:.3f}" for v in VIEWS) + " |\n")

    L.append("\n## (B) strong-vs-weak SD (per-organism ΔG terciles) — MCC per organism\n\n")
    L.append("| organism | n_val | ΔG cuts | " + " | ".join(VIEWS) + " |\n")
    L.append("|---|---|---|" + "|".join(["---"]*len(VIEWS)) + "|\n")
    for o in results["binary"]:
        r = results["binary"][o]
        L.append(f"| {o} | {r['n_val']} | {r['dG_tercile_cuts']} | " +
                 " | ".join(f"{r[v]:.3f}" for v in VIEWS) + " |\n")

    with open(os.path.join(args.out, "phase0_sd_report.md"), "w") as fh:
        fh.write("".join(L))
    print("".join(L))
    print(f"\nwrote {args.out}/phase0_sd_report.md + phase0_sd_results.json")


if __name__ == "__main__":
    main()
