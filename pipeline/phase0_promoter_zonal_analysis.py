#!/usr/bin/env python3
"""
Analyze promoter ZONAL features (robustness test #1) — is the promoter-vs-intergenic
signal localized to the core-promoter region, or carried by a co-occurring feature
elsewhere in the window?

For each zone, train a logistic-regression probe (promoter vs intergenic decoy) on the
committed train split using ONLY that zone's pooled features, and score MCC on val.

Views per zone: routing_<z> (96), emb_<z> (768), concat_<z> (864), kmer4_<z> (256,
from the zone's DNA), gc_<z> (1). 'full' reproduces the original whole-window pool.

Decision:
  localized   = MCC(core) − max(MCC(down), MCC(farup)) ≥ MARGIN  → signal is AT the
                promoter region; a mystery element elsewhere is NOT what drives it.
  p1_at_core  = routing_core > emb_core                          → MoE channel carries
                the promoter motif beyond the dense-accessible stream.

Output: <out>/phase0_promoter_zonal_report.md + phase0_promoter_zonal_results.json
Run: .venv/bin/python pipeline/phase0_promoter_zonal_analysis.py \
        --zonal data/datasets/phase0/phase0_features_promoters_zonal.npz \
        --all data/datasets/promoters/ALL.tsv --out experiments/phase0
"""
import argparse, csv, itertools, json, os
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import matthews_corrcoef

ZONE_BP = {"full": None, "minus35": (112, 124), "minus10": (136, 146),
           "core": (110, 150), "tss": (146, 156), "down": (156, 300), "farup": (0, 110)}
ZONE_ORDER = ["full", "core", "minus35", "minus10", "tss", "down", "farup"]
MARGIN = 0.10
_KIDX = {"".join(p): i for i, p in enumerate(itertools.product("ACGT", repeat=4))}


def kmer4(seq):
    v = np.zeros(256, np.float32); n = len(seq) - 3
    if n <= 0: return v
    c = 0
    for i in range(n):
        j = _KIDX.get(seq[i:i+4])
        if j is not None: v[j] += 1; c += 1
    return v / c if c else v


def gc(seq):
    s = seq.upper(); return (s.count("G")+s.count("C"))/len(s) if s else 0.0


def fit_mcc(X, y, tr, te):
    sc = StandardScaler(); Xtr = sc.fit_transform(X[tr]); Xte = sc.transform(X[te])
    clf = LogisticRegression(max_iter=5000, class_weight="balanced", n_jobs=-1, C=1.0)
    clf.fit(Xtr, y[tr])
    return round(float(matthews_corrcoef(y[te], clf.predict(Xte))), 4)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zonal", required=True)
    ap.add_argument("--all", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    d = np.load(args.zonal, allow_pickle=False)
    Z = {k: d[k] for k in d.files}
    ids = Z["ids"].astype(str); split = Z["split"].astype(str)
    y = (Z["label"].astype(str) == "positive").astype(int)
    tr = split == "train"; te = split == "val"

    seqmap = {}
    with open(args.all) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            seqmap[r["id"]] = r["window_seq"]
    seqs = [seqmap.get(i, "") for i in ids]

    rows = {}
    for z in ZONE_ORDER:
        emb, rout = Z[f"emb_{z}"], Z[f"routing_{z}"]
        concat = np.concatenate([emb, rout], axis=1)
        bp = ZONE_BP[z]
        zs = [s if bp is None else s[bp[0]:bp[1]] for s in seqs]
        km = np.stack([kmer4(s.upper()) for s in zs]).astype(np.float32)
        gcv = np.array([gc(s) for s in zs], np.float32).reshape(-1, 1)
        rows[z] = {
            "gc": fit_mcc(gcv, y, tr, te),
            "kmer4": fit_mcc(km, y, tr, te),
            "embedding": fit_mcc(emb, y, tr, te),
            "routing": fit_mcc(rout, y, tr, te),
            "concat": fit_mcc(concat, y, tr, te),
        }

    core = rows["core"]["concat"]
    loc_margin = round(core - max(rows["down"]["concat"], rows["farup"]["concat"]), 4)
    localized = loc_margin >= MARGIN
    p1_core = rows["core"]["routing"] > rows["core"]["embedding"]

    results = {"per_zone": rows,
               "decision": {"localized_to_core_promoter": bool(localized),
                            "localization_margin_mcc": loc_margin,
                            "p1_routing_gt_embedding_at_core": bool(p1_core),
                            "n_val": int(te.sum()), "n_val_pos": int(y[te].sum())}}
    with open(os.path.join(args.out, "phase0_promoter_zonal_results.json"), "w") as fh:
        json.dump(results, fh, indent=2)

    views = ["gc", "kmer4", "embedding", "routing", "concat"]
    L = ["# Phase-0 promoter robustness #1 — zonal (per-position) localization\n",
         "\nLogistic probe (promoter vs intergenic decoy) trained per ZONE on the committed\n"
         "train split; val MCC. Localizes WHERE in the 300-bp window the discriminative\n"
         "signal lives — to rule out a co-occurring confound outside the promoter.\n",
         f"\nval n={results['decision']['n_val']} ({results['decision']['n_val_pos']} promoters)\n",
         "\n## val MCC by zone × feature view\n\n",
         "| zone (bp) | " + " | ".join(views) + " |\n",
         "|---|" + "|".join(["---"]*len(views)) + "|\n"]
    for z in ZONE_ORDER:
        bp = ZONE_BP[z]; tag = "all" if bp is None else f"{bp[0]}–{bp[1]}"
        L.append(f"| {z} ({tag}) | " + " | ".join(f"{rows[z][v]:.3f}" for v in views) + " |\n")
    L.append(f"\n## Decision\n\n")
    L.append(f"- **Localized to core promoter: {'YES' if localized else 'NO'}** "
             f"(concat MCC core {core:.3f} − max(down {rows['down']['concat']:.3f}, "
             f"farup {rows['farup']['concat']:.3f}) = {loc_margin:+.3f}; margin ≥ {MARGIN})\n")
    L.append(f"- P1 at core (routing > embedding): {'YES' if p1_core else 'NO'} "
             f"({rows['core']['routing']:.3f} vs {rows['core']['embedding']:.3f})\n")
    L.append("\nReading: if `core` ≫ `down`/`farup`, the promoter call is driven by the\n"
             "−35/−10/TSS region, not by a co-occurring element elsewhere in the window.\n")
    with open(os.path.join(args.out, "phase0_promoter_zonal_report.md"), "w") as fh:
        fh.write("".join(L))
    print("".join(L))
    print(f"wrote {args.out}/phase0_promoter_zonal_report.md + .json")


if __name__ == "__main__":
    main()
