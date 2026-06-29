#!/usr/bin/env python3
"""
Phase-0 feature extraction (research_plan §4, §7) — runs ON the LBL cluster inside
the NeMo container, on the frozen GenomeOcean-MoE (pilot4 lr_match_dropout05).

For every 300-bp window in an element's ALL.tsv, do ONE frozen forward pass and
capture both feature views the Phase-0 smoke test needs:
  - embedding_only : 768-d mean-pooled last-layer hidden state (a dense model also
                     has this — the P1 control)
  - routing        : 96-d mean router softmax over 12 layers x 8 experts (the MoE
                     channel a dense model CANNOT supply)
`routing_concat` (864-d) is just np.concatenate([embedding, routing]) downstream.

Also stored (no model needed, computed here for convenience so the analysis stage is
self-contained and reproducible):
  - gc    : per-window G+C fraction      -> the GC-content-matched chance baseline (P0)
  - kmer4 : 256-d normalized 4-mer freqs -> a stronger pure-composition baseline

Reuses Junho's extract_features_moe + load_moe_model verbatim (no model code changes).

Output: <output-dir>/phase0_features_<element>.npz with arrays:
  emb (N,768) float32 | routing (N,96) float32 | gc (N,) float32 | kmer4 (N,256) float32
  ids, organism, phylum, label, split (all (N,) str)

Run (inside container, see submit_phase0.sh):
  python phase0_extract_features.py --element rbs \
     --checkpoint <ckpt> --config configs/config_100m_moe.yaml \
     --data-dir data/phase0 --output-dir experiments/results/phase0
"""
import argparse, csv, os, sys
import numpy as np
import torch

# make the eval/ sibling modules importable (same convention as eval_propose_fold_panel.py)
EVAL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "eval")
sys.path.insert(0, EVAL_DIR)
from eval_downstream import extract_features_moe          # noqa: E402
from eval_utils import (init_distributed, init_parallel_state,  # noqa: E402
                        load_moe_model)

ELEMENTS = ("promoters", "rbs", "rho")
_K = 4
_BASES = "ACGT"
_KMER_INDEX = {"".join(p): i for i, p in enumerate(
    __import__("itertools").product(_BASES, repeat=_K))}


def gc_fraction(seq):
    s = seq.upper()
    n = len(s)
    if n == 0:
        return 0.0
    return (s.count("G") + s.count("C")) / n


def kmer4_vector(seq):
    s = seq.upper()
    v = np.zeros(len(_KMER_INDEX), dtype=np.float32)
    n = len(s) - _K + 1
    if n <= 0:
        return v
    cnt = 0
    for i in range(n):
        idx = _KMER_INDEX.get(s[i:i + _K])
        if idx is not None:        # skip windows of N etc.
            v[idx] += 1.0
            cnt += 1
    if cnt:
        v /= cnt
    return v


def load_windows(data_dir, element):
    """Read ALL.tsv + split_80_20.tsv for an element; return joined list of dicts."""
    all_tsv   = os.path.join(data_dir, f"{element}_ALL.tsv")
    split_tsv = os.path.join(data_dir, f"{element}_split_80_20.tsv")

    split_of = {}
    with open(split_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            split_of[r["window_id"]] = r["split"]

    rows = []
    with open(all_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            wid = r["id"]
            rows.append(dict(
                id=wid,
                organism=r["organism"],
                phylum=r["phylum"],
                label=r["label"],
                split=split_of.get(wid, "none"),
                seq=r["window_seq"],
            ))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--element", required=True, choices=ELEMENTS + ("all",))
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--tokenizer", default="DOEJGI/GenomeOcean-4B")
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--seq-length", type=int, default=1024)   # model context (RoPE)
    ap.add_argument("--max-length", type=int, default=128)    # tokens per 300-bp window (~61 used)
    ap.add_argument("--batch-size", type=int, default=64)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    elements = ELEMENTS if args.element == "all" else (args.element,)

    print("Initializing Megatron parallel state (TP=PP=EP=1) ...", flush=True)
    init_distributed()
    init_parallel_state()

    print(f"Loading frozen MoE from {args.checkpoint} ...", flush=True)
    model, tokenizer = load_moe_model(args)
    device = next(model.parameters()).device
    print(f"  model loaded on {device}", flush=True)

    for element in elements:
        rows = load_windows(args.data_dir, element)
        seqs = [r["seq"] for r in rows]
        print(f"\n[{element}] {len(seqs):,} windows — extracting features ...", flush=True)

        emb, routing, _dom = extract_features_moe(
            model, tokenizer, seqs,
            batch_size=args.batch_size, max_length=args.max_length, device=device)

        print(f"[{element}] emb={emb.shape} routing={routing.shape} — "
              f"computing GC + 4-mer baselines ...", flush=True)
        gc    = np.array([gc_fraction(s) for s in seqs], dtype=np.float32)
        kmer4 = np.stack([kmer4_vector(s) for s in seqs]).astype(np.float32)

        out = os.path.join(args.output_dir, f"phase0_features_{element}.npz")
        np.savez_compressed(
            out,
            emb=emb.astype(np.float32),
            routing=routing.astype(np.float32),
            gc=gc, kmer4=kmer4,
            ids=np.array([r["id"]       for r in rows]),
            organism=np.array([r["organism"] for r in rows]),
            phylum=np.array([r["phylum"]   for r in rows]),
            label=np.array([r["label"]    for r in rows]),
            split=np.array([r["split"]    for r in rows]),
        )
        print(f"[{element}] wrote {out}", flush=True)


if __name__ == "__main__":
    main()
