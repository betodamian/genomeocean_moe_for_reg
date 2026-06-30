#!/usr/bin/env python3
"""
Extract frozen-MoE features for the confound-free promoter-detection test
(within-class dinucleotide-shuffle). Runs ON the LBL cluster.

Reads data/phase0/promoter_dinuc_ALL.tsv (id, organism, split, klass, seq;
klass in {P, P_shuf, I, I_shuf}) and does one frozen forward pass per row:
  emb (768)  + routing (96)  per window.

Downstream analysis (phase0_promoter_dinuc_analysis.py) computes
  MCC(P vs P_shuf) − MCC(I vs I_shuf)  per feature view (the confound-free signal).

Output: <out>/phase0_features_promoter_dinuc.npz
Run inside container — see submit_phase0_promoter_dinuc.sh.
"""
import argparse, csv, os, sys
import numpy as np

EVAL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "eval")
sys.path.insert(0, EVAL_DIR)
from eval_downstream import extract_features_moe                    # noqa: E402
from eval_utils import init_distributed, init_parallel_state, load_moe_model  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--tokenizer", default="DOEJGI/GenomeOcean-4B")
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--seq-length", type=int, default=1024)
    ap.add_argument("--max-length", type=int, default=64)   # 81 bp ~ 17 tokens
    ap.add_argument("--batch-size", type=int, default=64)
    args = ap.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    rows = []
    with open(os.path.join(args.data_dir, "promoter_dinuc_ALL.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)
    seqs = [r["seq"] for r in rows]
    print(f"promoter dinuc-test windows: {len(seqs):,}", flush=True)

    init_distributed(); init_parallel_state()
    print(f"loading frozen MoE from {args.checkpoint} ...", flush=True)
    model, tokenizer = load_moe_model(args)
    device = next(model.parameters()).device
    print(f"  model on {device}; extracting ...", flush=True)

    emb, routing, _dom = extract_features_moe(
        model, tokenizer, seqs, batch_size=args.batch_size,
        max_length=args.max_length, device=device)

    out = os.path.join(args.output_dir, "phase0_features_promoter_dinuc.npz")
    np.savez_compressed(
        out,
        emb=emb.astype(np.float32), routing=routing.astype(np.float32),
        ids=np.array([r["id"] for r in rows]),
        organism=np.array([r["organism"] for r in rows]),
        split=np.array([r["split"] for r in rows]),
        klass=np.array([r["klass"] for r in rows]))
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
