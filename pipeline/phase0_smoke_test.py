#!/usr/bin/env python3
"""
Phase-0 smoke test / ceiling benchmark (research_plan §4, P0) — the go/no-go gate.

Question (P0): on SEEN data, does the frozen model's `routing_concat` discriminate
each element from its SAME-CONTEXT decoy ABOVE a GC-content-matched chance baseline?

Reads the NPZ feature files from phase0_extract_features.py and, for each
pre-registered task, trains a logistic-regression probe on the committed 80/20
cluster-level split and evaluates on validation.

Feature views (P1 controls):
  embedding_only (768)  — a dense model also has this
  routing_only   (96)   — the pure MoE routing channel
  routing_concat (864)  — PRIMARY detector (embedding + routing)
Baselines:
  gc_only (1)   — GC fraction; the pre-registered chance baseline (P0)
  kmer4   (256) — 4-mer composition; a stronger "no-model" sequence baseline

Tasks (same-context decoys, organism-confound controlled):
  promoter_vs_intergenic    promoter        vs intergenic (per-organism matched)
  rbs_TIS_vs_intergenic     SD+UNSD starts  vs intergenic (per-organism matched)
  rbs_SD_vs_UNSD            SD starts       vs UNSD starts (Tier-2 PRIMARY; GC-clean)
  rho_t1_vs_intergenic      in-vivo Rho     vs intergenic (per-organism matched)
  rho_t1_vs_intrinsic_ecoli Rho (E.coli)    vs intrinsic terminators (Tier-2, E.coli-only)

Metrics: MCC (primary), AUROC, AUPRC, F1, accuracy + bootstrap 95% CI on MCC & AUPRC.
GATE per task: routing_concat passes if its MCC lower-CI > gc_only MCC upper-CI
(i.e. it clears the GC baseline with non-overlapping 95% CIs).

Output: <out>/phase0_report.md + phase0_results.json
Run:    .venv/bin/python pipeline/phase0_smoke_test.py --features-dir <dir> --out <dir>
"""
import argparse, json, os
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (matthews_corrcoef, f1_score, accuracy_score,
                             roc_auc_score, average_precision_score)

RNG = np.random.default_rng(20260629)
N_BOOT = 1000

# task -> (element, positive_labels, negative_labels, restrict_organism_or_None)
TASKS = {
    "promoter_vs_intergenic": ("promoters",
        {"positive"}, {"decoy"}, None),
    "rbs_TIS_vs_intergenic": ("rbs",
        {"positive_SD", "positive_UNSD"}, {"decoy_intergenic"}, None),
    "rbs_SD_vs_UNSD": ("rbs",
        {"positive_SD"}, {"positive_UNSD"}, None),
    "rho_t1_vs_intergenic": ("rho",
        {"positive_rho_t1", "positive_rho_t1_rsr"}, {"decoy_intergenic"}, None),
    "rho_t1_vs_intrinsic_ecoli": ("rho",
        {"positive_rho_t1", "positive_rho_t1_rsr"}, {"decoy_intrinsic"},
        "ecoli_K12_MG1655"),
}

FEATURE_VIEWS = ["gc_only", "kmer4", "embedding_only", "routing_only", "routing_concat"]


def load_npz(features_dir, element):
    p = os.path.join(features_dir, f"phase0_features_{element}.npz")
    d = np.load(p, allow_pickle=False)
    return {k: d[k] for k in d.files}


def build_view(data, view):
    if view == "gc_only":        return data["gc"].reshape(-1, 1)
    if view == "kmer4":          return data["kmer4"]
    if view == "embedding_only": return data["emb"]
    if view == "routing_only":   return data["routing"]
    if view == "routing_concat": return np.concatenate([data["emb"], data["routing"]], axis=1)
    raise ValueError(view)


def bootstrap_ci(y_true, y_pred, y_score, metric, n=N_BOOT):
    """95% bootstrap CI for a metric over the validation set."""
    N = len(y_true)
    vals = []
    for _ in range(n):
        idx = RNG.integers(0, N, N)
        yt = y_true[idx]
        if len(np.unique(yt)) < 2:
            continue
        if metric == "mcc":
            vals.append(matthews_corrcoef(yt, y_pred[idx]))
        elif metric == "auprc":
            vals.append(average_precision_score(yt, y_score[idx]))
    if not vals:
        return (float("nan"), float("nan"))
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def evaluate(Xtr, ytr, Xte, yte, view):
    scaler = StandardScaler()
    Xtr_s = scaler.fit_transform(Xtr)
    Xte_s = scaler.transform(Xte)
    clf = LogisticRegression(max_iter=5000, solver="lbfgs",
                             class_weight="balanced", n_jobs=-1, C=1.0)
    clf.fit(Xtr_s, ytr)
    y_pred  = clf.predict(Xte_s)
    y_score = clf.predict_proba(Xte_s)[:, 1]

    mcc_lo, mcc_hi = bootstrap_ci(yte, y_pred, y_score, "mcc")
    ap_lo,  ap_hi  = bootstrap_ci(yte, y_pred, y_score, "auprc")
    return {
        "mcc": round(matthews_corrcoef(yte, y_pred), 4),
        "mcc_ci": [round(mcc_lo, 4), round(mcc_hi, 4)],
        "auroc": round(roc_auc_score(yte, y_score), 4),
        "auprc": round(average_precision_score(yte, y_score), 4),
        "auprc_ci": [round(ap_lo, 4), round(ap_hi, 4)],
        "f1": round(f1_score(yte, y_pred), 4),
        "accuracy": round(accuracy_score(yte, y_pred), 4),
    }


def per_organism_mcc(data, keep, split, pos_labels, ytr, tr_mask):
    """routing_concat val MCC broken down by organism (confound visibility).

    One probe trained on all task training data, then evaluated within each organism
    on the validation set — so we can see if a pooled win is organism-driven.
    """
    Xall = build_view(data, "routing_concat")
    scaler = StandardScaler(); Xtr_s = scaler.fit_transform(Xall[tr_mask])
    clf = LogisticRegression(max_iter=5000, solver="lbfgs",
                             class_weight="balanced", n_jobs=-1, C=1.0)
    clf.fit(Xtr_s, ytr)
    out = {}
    val_mask = keep & (split == "val")
    for org in sorted(set(data["organism"][val_mask])):
        om = val_mask & (data["organism"] == org)
        if om.sum() < 10:
            continue
        yt = np.isin(data["label"][om], list(pos_labels)).astype(int)
        if len(np.unique(yt)) < 2:
            continue
        yp = clf.predict(scaler.transform(Xall[om]))
        out[org] = {"n": int(om.sum()), "mcc": round(matthews_corrcoef(yt, yp), 4)}
    return out


def run_task(task, cfg, features_dir):
    element, pos_labels, neg_labels, restrict_org = cfg
    data = load_npz(features_dir, element)

    keep = np.isin(data["label"], list(pos_labels | neg_labels))
    if restrict_org:
        keep &= (data["organism"] == restrict_org)

    y = np.isin(data["label"], list(pos_labels)).astype(int)
    split = data["split"]
    tr_mask = keep & (split == "train")
    te_mask = keep & (split == "val")

    info = {
        "element": element,
        "positives": sorted(pos_labels), "negatives": sorted(neg_labels),
        "restrict_organism": restrict_org,
        "n_train": int(tr_mask.sum()), "n_val": int(te_mask.sum()),
        "n_train_pos": int(y[tr_mask].sum()), "n_val_pos": int(y[te_mask].sum()),
        "views": {},
    }
    if info["n_train_pos"] < 5 or info["n_val_pos"] < 5 or \
       (info["n_train"] - info["n_train_pos"]) < 5 or (info["n_val"] - info["n_val_pos"]) < 5:
        info["error"] = "insufficient class counts"
        return info

    ytr, yte = y[tr_mask], y[te_mask]
    for view in FEATURE_VIEWS:
        X = build_view(data, view)
        info["views"][view] = evaluate(X[tr_mask], ytr, X[te_mask], yte, view)

    # gate: routing_concat MCC lower-CI > gc_only MCC upper-CI
    rc = info["views"]["routing_concat"]
    gc = info["views"]["gc_only"]
    info["gate_pass"] = bool(rc["mcc_ci"][0] > gc["mcc_ci"][1])
    info["gate_margin_mcc"] = round(rc["mcc"] - gc["mcc"], 4)

    # per-organism breakdown for routing_concat (confound visibility)
    info["routing_concat_per_organism"] = per_organism_mcc(
        data, keep, split, pos_labels, ytr, tr_mask)
    return info


def fmt_view(v):
    return (f"MCC={v['mcc']:.3f} [{v['mcc_ci'][0]:.3f},{v['mcc_ci'][1]:.3f}]  "
            f"AUPRC={v['auprc']:.3f}  AUROC={v['auroc']:.3f}  F1={v['f1']:.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--features-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    results = {}
    for task, cfg in TASKS.items():
        print(f"running {task} ...", flush=True)
        results[task] = run_task(task, cfg, args.features_dir)

    with open(os.path.join(args.out, "phase0_results.json"), "w") as fh:
        json.dump(results, fh, indent=2)

    # ── markdown report ────────────────────────────────────────────────────────
    L = ["# Phase-0 smoke test — ceiling benchmark (research_plan §4, P0)\n",
         "Frozen GenomeOcean-MoE (pilot4 lr_match_dropout05). Linear probe on the "
         "committed 80/20 cluster-level split. **Gate:** `routing_concat` MCC "
         "lower-CI > `gc_only` MCC upper-CI (clears the GC-matched baseline with "
         "non-overlapping 95% CIs).\n"]

    n_pass = sum(1 for r in results.values() if r.get("gate_pass"))
    L.append(f"\n## Verdict: {n_pass}/{len(results)} tasks PASS the GC-matched gate\n")
    L.append("\n| task | n_val (pos) | routing_concat MCC | gc_only MCC | gate |\n")
    L.append("|---|---|---|---|---|\n")
    for t, r in results.items():
        if "error" in r:
            L.append(f"| {t} | — | — | — | ERROR: {r['error']} |\n"); continue
        rc, gc = r["views"]["routing_concat"], r["views"]["gc_only"]
        gate = "**PASS**" if r["gate_pass"] else "FAIL"
        L.append(f"| {t} | {r['n_val']} ({r['n_val_pos']}) | "
                 f"{rc['mcc']:.3f} [{rc['mcc_ci'][0]:.3f},{rc['mcc_ci'][1]:.3f}] | "
                 f"{gc['mcc']:.3f} [{gc['mcc_ci'][0]:.3f},{gc['mcc_ci'][1]:.3f}] | {gate} |\n")

    for t, r in results.items():
        L.append(f"\n## {t}\n")
        if "error" in r:
            L.append(f"\n**ERROR:** {r['error']}\n"); continue
        L.append(f"\n- element: `{r['element']}`  positives={r['positives']}  "
                 f"negatives={r['negatives']}  restrict_org={r['restrict_organism']}\n")
        L.append(f"- train: {r['n_train']} ({r['n_train_pos']} pos)  "
                 f"val: {r['n_val']} ({r['n_val_pos']} pos)\n")
        L.append(f"- gate: {'PASS' if r['gate_pass'] else 'FAIL'}  "
                 f"(routing_concat − gc_only MCC = {r['gate_margin_mcc']:+.3f})\n\n")
        L.append("| feature view | MCC [95% CI] | AUPRC | AUROC | F1 | acc |\n")
        L.append("|---|---|---|---|---|---|\n")
        for view in FEATURE_VIEWS:
            v = r["views"][view]
            L.append(f"| {view} | {v['mcc']:.3f} [{v['mcc_ci'][0]:.3f},{v['mcc_ci'][1]:.3f}] | "
                     f"{v['auprc']:.3f} | {v['auroc']:.3f} | {v['f1']:.3f} | {v['accuracy']:.3f} |\n")
        po = r.get("routing_concat_per_organism", {})
        if po:
            L.append("\nrouting_concat val MCC by organism: " +
                     ", ".join(f"{o}={d['mcc']:.3f}(n={d['n']})" for o, d in po.items()) + "\n")

    with open(os.path.join(args.out, "phase0_report.md"), "w") as fh:
        fh.write("".join(L))
    print("\n" + "".join(L))
    print(f"\nwrote {args.out}/phase0_report.md + phase0_results.json")


if __name__ == "__main__":
    main()
