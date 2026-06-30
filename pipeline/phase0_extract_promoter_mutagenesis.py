#!/usr/bin/env python3
"""
Phase-0 promoter ROBUSTNESS test #2 — in-silico motif mutagenesis (causal control,
research_plan §4 / §9b). Runs ON the LBL cluster, frozen MoE.

The strongest test that the model uses the PROMOTER and not a co-occurring "mystery"
element: take promoter windows the model detects, GC-PRESERVINGLY scramble the
promoter motif, and re-run the frozen forward pass. If destroying the motif collapses
the model's promoter prediction while scrambling a MATCHED control region (same size,
downstream, no promoter) does not, the model is causally using the promoter.

Window geometry: 300 bp, TSS at index 150 (build_promoter_windows.py). Variants
(GC-preserving in-place shuffle = exact base composition kept, so any prediction drop
is NOT a GC artifact; per-window deterministic RNG):
  original   — unperturbed (sanity: should score like the original probe)
  m35        — scramble −35 hexamer  [115,121)
  m10        — scramble −10 hexamer  [138,144)
  m10m35     — scramble both hexamers
  core       — scramble core promoter [110,148)  (−40..−2; covers archaeal TATA/BRE too)
  ctrl       — scramble matched downstream [160,198)  (+10..+48, same 38 bp) → CONTROL

Only POSITIVE windows are perturbed (we test whether detection survives motif loss).
The probe itself is trained downstream on the ORIGINAL features (promoter vs decoy).

Decision (phase0_promoter_mutagenesis_analysis.py):
  ΔP(promoter) for core / m10m35  ≫  ΔP for ctrl  → causal use of the promoter motif.
  (bacteria should respond to m10/m35; archaeal H. volcanii to core not bacterial
   hexamers — a built-in specificity check.)

Output: <out>/phase0_features_promoters_mutagenesis.npz  (emb_<v>, routing_<v> per variant)
Run inside container — see submit_phase0_promoter_robustness.sh.
"""
import argparse, csv, os, sys
import numpy as np

EVAL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "eval")
sys.path.insert(0, EVAL_DIR)
from eval_downstream import extract_features_moe                    # noqa: E402
from eval_utils import init_distributed, init_parallel_state, load_moe_model  # noqa: E402

# variant -> list of (lo, hi) regions to GC-preservingly scramble (window indices)
VARIANTS = {
    "original": [],
    "m35":      [(115, 121)],
    "m10":      [(138, 144)],
    "m10m35":   [(115, 121), (138, 144)],
    "core":     [(110, 148)],
    "ctrl":     [(160, 198)],
}


def scramble(seq, regions, seed):
    if not regions:
        return seq
    rng = np.random.default_rng(seed)
    s = list(seq)
    for lo, hi in regions:
        lo, hi = max(0, lo), min(len(s), hi)
        seg = s[lo:hi]
        perm = rng.permutation(len(seg))
        s[lo:hi] = [seg[i] for i in perm]   # exact composition preserved
    return "".join(s)


def load_promoter_positives(data_dir, split_only=None):
    all_tsv   = os.path.join(data_dir, "promoters_ALL.tsv")
    split_tsv = os.path.join(data_dir, "promoters_split_80_20.tsv")
    split_of = {}
    with open(split_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            split_of[r["window_id"]] = r["split"]
    rows = []
    with open(all_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["label"] != "positive":
                continue
            sp = split_of.get(r["id"], "none")
            if split_only and sp != split_only:
                continue
            rows.append(dict(id=r["id"], organism=r["organism"], phylum=r["phylum"],
                             split=sp, seq=r["window_seq"]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--tokenizer", default="DOEJGI/GenomeOcean-4B")
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--split", default="val", help="which split's positives to perturb (cheap: val)")
    ap.add_argument("--seq-length", type=int, default=1024)
    ap.add_argument("--max-length", type=int, default=128)
    ap.add_argument("--batch-size", type=int, default=64)
    args = ap.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    rows = load_promoter_positives(args.data_dir, split_only=args.split)
    print(f"promoter positives ({args.split}): {len(rows):,}", flush=True)

    init_distributed(); init_parallel_state()
    print(f"loading frozen MoE from {args.checkpoint} ...", flush=True)
    model, tokenizer = load_moe_model(args)
    device = next(model.parameters()).device

    save = {"ids": np.array([r["id"] for r in rows]),
            "organism": np.array([r["organism"] for r in rows]),
            "phylum": np.array([r["phylum"] for r in rows]),
            "split": np.array([r["split"] for r in rows])}

    for v, regions in VARIANTS.items():
        seqs = [scramble(r["seq"], regions, seed=(i * 7919 + 13))
                for i, r in enumerate(rows)]
        print(f"[{v}] extracting {len(seqs):,} windows ...", flush=True)
        emb, routing, _dom = extract_features_moe(
            model, tokenizer, seqs, batch_size=args.batch_size,
            max_length=args.max_length, device=device)
        save[f"emb_{v}"] = emb.astype(np.float32)
        save[f"routing_{v}"] = routing.astype(np.float32)

    out = os.path.join(args.output_dir, "phase0_features_promoters_mutagenesis.npz")
    np.savez_compressed(out, **save)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
