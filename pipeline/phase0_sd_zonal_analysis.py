#!/usr/bin/env python3
"""
Analyze the zonal (per-position) SD features — settle pooling-dilution (a) vs
non-specialization (b).

Within-organism ΔG regression (Spearman ρ) across feature views:
  gc_sd        — GC of the SD region (1-d)            no-model baseline
  kmer4_sd     — 4-mer freqs of the SD region (256-d) no-model composition baseline
  emb_full     — embedding, whole-window mean-pool    (the original view)
  routing_full — routing,    whole-window mean-pool    (original; ρ≈0.05 in pooled run)
  emb_sd       — embedding pooled over SD-zone tokens ONLY
  routing_sd   — routing   pooled over SD-zone tokens ONLY   <- the decisive view
  routing_start— routing   pooled over the start-codon token
  concat_sd    — emb_sd + routing_sd

Decision rule:
  routing_sd >> routing_full        -> (a) dilution: signal was there, un-pooling found it
  routing_sd ~ 0 (≈ routing_full)   -> (b) non-specialization: router doesn't encode SD
  (and P1 at the SD zone: routing_sd > emb_sd ?)

Output: <out>/phase0_sd_zonal_report.md + phase0_sd_zonal_results.json
Run: .venv/bin/python pipeline/phase0_sd_zonal_analysis.py \
        --zonal data/datasets/phase0/phase0_features_rbs_zonal.npz \
        --deltaG data/datasets/rbs/sd_deltaG.tsv \
        --all data/datasets/rbs/ALL.tsv --out experiments/phase0
"""
import argparse, csv, itertools, json, os
import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

SD0, SD1 = 128, 148        # SD-zone bp window (matches extractor)
MIN_N = 40
_KIDX = {"".join(p): i for i, p in enumerate(itertools.product("ACGT", repeat=4))}


def kmer4(seq):
    v = np.zeros(256, np.float32); n = len(seq) - 3
    if n <= 0: return v
    c = 0
    for i in range(n):
        j = _KIDX.get(seq[i:i+4])
        if j is not None: v[j] += 1; c += 1
    return v / c if c else v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zonal", required=True)
    ap.add_argument("--deltaG", required=True)
    ap.add_argument("--all", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    d = np.load(args.zonal, allow_pickle=False)
    Z = {k: d[k] for k in d.files}
    ids = Z["ids"].astype(str); org = Z["organism"].astype(str); split = Z["split"].astype(str)

    dgmap = {}
    with open(args.deltaG) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            dgmap[r["window_id"]] = float(r["dG"])

    seqmap = {}
    with open(args.all) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["label"].startswith("positive_"):
                seqmap[r["id"]] = r["window_seq"]

    dG = np.array([dgmap.get(i, np.nan) for i in ids])
    sd_seq = [seqmap.get(i, "")[SD0:SD1].upper() for i in ids]
    gc_sd = np.array([( (s.count("G")+s.count("C"))/len(s) if s else 0.0) for s in sd_seq], np.float32)
    km_sd = np.stack([kmer4(s) for s in sd_seq]).astype(np.float32)

    views = {
        "gc_sd": gc_sd.reshape(-1, 1),
        "kmer4_sd": km_sd,
        "emb_full": Z["emb_full"], "routing_full": Z["routing_full"],
        "emb_sd": Z["emb_sd"], "routing_sd": Z["routing_sd"],
        "routing_start": Z["routing_start"],
        "concat_sd": np.concatenate([Z["emb_sd"], Z["routing_sd"]], axis=1),
    }
    VIEW_ORDER = list(views)

    has = ~np.isnan(dG)
    organisms = sorted(set(org[has]))
    per_org = {}
    for o in organisms:
        om = has & (org == o)
        tr = om & (split == "train"); te = om & (split == "val")
        if tr.sum() < MIN_N or te.sum() < MIN_N:
            continue
        res = {"n_val": int(te.sum())}
        for v in VIEW_ORDER:
            X = views[v]
            sc = StandardScaler(); Xtr = sc.fit_transform(X[tr]); Xte = sc.transform(X[te])
            rho = spearmanr(dG[te], Ridge(alpha=1.0).fit(Xtr, dG[tr]).predict(Xte)).correlation
            res[v] = round(float(rho), 4)
        per_org[o] = res

    macro = {v: round(float(np.mean([per_org[o][v] for o in per_org])), 4) for v in VIEW_ORDER}

    # decision
    dilution = macro["routing_sd"] > macro["routing_full"] + 0.05
    p1_sd = macro["routing_sd"] > macro["emb_sd"]
    beats_kmer = macro["concat_sd"] > macro["kmer4_sd"]
    branch = ("(a) pooling dilution — SD-zone routing recovers signal lost to whole-window pooling"
              if dilution else
              "(b) non-specialization — un-pooling does NOT recover SD; router does not encode it")

    results = {"per_organism": per_org, "macro": macro,
               "decision": {"branch": branch, "dilution": bool(dilution),
                            "p1_at_sd_zone": bool(p1_sd), "concat_sd_beats_kmer4_sd": bool(beats_kmer)}}
    with open(os.path.join(args.out, "phase0_sd_zonal_results.json"), "w") as fh:
        json.dump(results, fh, indent=2)

    L = ["# Phase-0 RBS SD — zonal (per-position) re-analysis\n",
         "\nWithin-organism ΔG regression (Spearman ρ). Routing/embedding pooled over the\n"
         "SD-zone tokens ONLY (no whole-window averaging), vs the original full-window pool.\n"]
    L.append("\n## Macro-average across organisms (Spearman ρ predicting SD ΔG)\n\n")
    L.append("| view | ρ |\n|---|---|\n")
    for v in VIEW_ORDER:
        L.append(f"| {v} | {macro[v]:.3f} |\n")
    L.append(f"\n**Decision: {branch}**\n\n")
    L.append(f"- routing_sd {macro['routing_sd']:.3f} vs routing_full {macro['routing_full']:.3f} "
             f"→ dilution {'YES' if dilution else 'NO'}\n")
    L.append(f"- P1 at SD zone (routing_sd > emb_sd): {'YES' if p1_sd else 'NO'} "
             f"({macro['routing_sd']:.3f} vs {macro['emb_sd']:.3f})\n")
    L.append(f"- concat_sd > kmer4_sd (model beats SD-region composition): "
             f"{'YES' if beats_kmer else 'NO'} ({macro['concat_sd']:.3f} vs {macro['kmer4_sd']:.3f})\n")
    L.append("\n## Per-organism ρ\n\n")
    L.append("| organism | n_val | " + " | ".join(VIEW_ORDER) + " |\n")
    L.append("|---|---|" + "|".join(["---"]*len(VIEW_ORDER)) + "|\n")
    for o in per_org:
        r = per_org[o]
        L.append(f"| {o} | {r['n_val']} | " + " | ".join(f"{r[v]:.3f}" for v in VIEW_ORDER) + " |\n")

    with open(os.path.join(args.out, "phase0_sd_zonal_report.md"), "w") as fh:
        fh.write("".join(L))
    print("".join(L))
    print(f"\nwrote {args.out}/phase0_sd_zonal_report.md + phase0_sd_zonal_results.json")


if __name__ == "__main__":
    main()
