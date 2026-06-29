#!/usr/bin/env python3
"""
Phase-0 RBS SD method fix (research_plan §4 mode-2/-3, §5d pre-registered sweep) —
the DECISIVE experiment to distinguish the two remaining SD hypotheses:

  (a) pooling dilution  — SD signal is real but mean-pooling over ~60 tokens washes
                          out the 1-2 token SD motif (PHASE0_FINDINGS.md addendum)
  (b) non-specialization — the router simply does not encode a short regulatory motif

Instead of mean-pooling routing/embedding over the WHOLE 300-bp window, this captures
PER-TOKEN routing + hidden states and pools them within ZONES defined relative to the
start codon (TIS at window index 150). Token→bp coordinates are reconstructed from BPE
token strings (each token is a literal DNA substring; method from
eval_expert_specialization.build_token_coordinate_map), so we pool over exactly the
tokens covering the SD region — no whole-window averaging.

Zones (bp in the 300-bp window; window[150] = first base of start codon):
  sd    = [128,148)  (-22..-3)  — where the Shine-Dalgarno sits
  start = [148,154)  (-2..+4)   — the start codon
  down  = [154,214)  (+4..+64)  — early ORF body
  full  = all content tokens     — reproduces the original mean-pool (alignment control)

If routing_sd >> routing_full → branch (a) dilution (signal was there, un-pooled it).
If routing_sd ~ 0 while emb_sd / kmer4 carry it → branch (b) non-specialization.

Output: <out>/phase0_features_rbs_zonal.npz with, for each zone z in {full,sd,start,down}:
  emb_<z> (N,768)  routing_<z> (N,96)   + ids/organism/phylum/label/split/gc
Run (inside container, see submit_phase0_zonal.sh).
"""
import argparse, csv, os, sys
import numpy as np
import torch

EVAL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "eval")
sys.path.insert(0, EVAL_DIR)
from eval_utils import init_distributed, init_parallel_state, load_moe_model  # noqa: E402

ZONES = {"sd": (128, 148), "start": (148, 154), "down": (154, 214)}  # full handled separately
TIS_IDX = 150


def load_rbs_positives(data_dir):
    all_tsv   = os.path.join(data_dir, "rbs_ALL.tsv")
    split_tsv = os.path.join(data_dir, "rbs_split_80_20.tsv")
    split_of = {}
    with open(split_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            split_of[r["window_id"]] = r["split"]
    rows = []
    with open(all_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r["label"].startswith("positive_"):
                continue
            rows.append(dict(id=r["id"], organism=r["organism"], phylum=r["phylum"],
                             label=r["label"], split=split_of.get(r["id"], "none"),
                             seq=r["window_seq"]))
    return rows


def token_bp_coords(token_strings):
    """(bp_start, bp_end) per token; each BPE token is a literal DNA substring."""
    coords, bp = [], 0
    for ts in token_strings:
        clean = ts.replace("▁", "").replace("Ġ", "").strip().upper()
        clean = "".join(c for c in clean if c in "ACGTN")
        coords.append((bp, bp + len(clean))); bp += len(clean)
    return coords


def zone_mask(coords, z0, z1):
    """boolean per-token: token overlaps [z0,z1)."""
    return np.array([(s < z1 and e > z0) for (s, e) in coords], dtype=bool)


def extract_zonal(model, tokenizer, sequences, batch_size, max_length, device):
    """Return dict zone -> (emb (N,768), routing (N,96)), zones = full + ZONES keys."""
    mc = model.module if hasattr(model, "module") else model
    cfg = mc.config if hasattr(mc, "config") else model.config
    num_layers, num_experts = cfg.num_layers, cfg.num_moe_experts

    captured = {}
    def hidden_hook(m, i, o): captured["hidden"] = o.detach()
    h_handle = mc.decoder.final_layernorm.register_forward_hook(hidden_hook)

    router_captured = {}
    def make_hook(idx):
        def hook(m, i, o):
            if isinstance(o, tuple) and len(o) >= 2 and isinstance(o[0], torch.Tensor) \
               and o[0].dim() == 2 and o[0].shape[-1] == num_experts:
                router_captured[idx] = o[0].float().detach()
        return hook
    r_handles = [layer.mlp.router.register_forward_hook(make_hook(i))
                 for i, layer in enumerate(mc.decoder.layers)
                 if hasattr(layer, "mlp") and hasattr(layer.mlp, "router")]

    special = {sid for sid in [getattr(tokenizer, "cls_token_id", None),
                               getattr(tokenizer, "sep_token_id", None),
                               getattr(tokenizer, "pad_token_id", None)] if sid is not None}

    zone_names = ["full"] + list(ZONES)
    out_emb = {z: [] for z in zone_names}
    out_rout = {z: [] for z in zone_names}
    n_mismatch = 0

    try:
        for start in range(0, len(sequences), batch_size):
            batch = sequences[start:start + batch_size]
            enc = tokenizer(batch, max_length=max_length, truncation=True,
                            padding="max_length", return_tensors="pt")
            input_ids = enc["input_ids"].to(device)
            pos = torch.arange(input_ids.shape[1], device=device).unsqueeze(0).expand_as(input_ids)
            router_captured.clear()
            with torch.no_grad():
                mc.forward(input_ids=input_ids, position_ids=pos, attention_mask=None)
            B, S = input_ids.shape
            hidden = captured["hidden"].transpose(0, 1)              # (B,S,768)
            # per-token routing: (B,S,layers*experts)
            parts = []
            for li in range(num_layers):
                if li in router_captured:
                    parts.append(router_captured[li].reshape(B, S, num_experts))
                else:
                    parts.append(torch.zeros(B, S, num_experts, device=device))
            routing = torch.cat(parts, dim=-1).cpu().float().numpy()  # (B,S,96)
            hidden_np = hidden.cpu().float().numpy()                  # (B,S,768)
            am = enc["attention_mask"].numpy().astype(bool)
            ids_np = input_ids.cpu().numpy()

            for b in range(B):
                content = am[b].copy()
                for sid in special:
                    content &= (ids_np[b] != sid)
                cidx = np.where(content)[0]
                toks = tokenizer.convert_ids_to_tokens(ids_np[b][cidx].tolist())
                coords = token_bp_coords(toks)
                emb_t = hidden_np[b][cidx]        # (T,768)
                rout_t = routing[b][cidx]         # (T,96)
                if len(coords) != emb_t.shape[0]:
                    n_mismatch += 1
                # full = mean over all content tokens
                out_emb["full"].append(emb_t.mean(0) if len(emb_t) else np.zeros(768, np.float32))
                out_rout["full"].append(rout_t.mean(0) if len(rout_t) else np.zeros(96, np.float32))
                for z, (z0, z1) in ZONES.items():
                    m = zone_mask(coords, z0, z1)
                    if m.any():
                        out_emb[z].append(emb_t[m].mean(0))
                        out_rout[z].append(rout_t[m].mean(0))
                    else:
                        out_emb[z].append(np.zeros(768, np.float32))
                        out_rout[z].append(np.zeros(96, np.float32))

            if (start // batch_size + 1) % 50 == 0:
                print(f"  {start+len(batch)}/{len(sequences)}", flush=True)
    finally:
        h_handle.remove()
        for h in r_handles: h.remove()

    print(f"  token/coord mismatches: {n_mismatch}/{len(sequences)}", flush=True)
    return ({z: np.stack(out_emb[z]).astype(np.float32) for z in zone_names},
            {z: np.stack(out_rout[z]).astype(np.float32) for z in zone_names})


def gc_fraction(s):
    s = s.upper(); n = len(s)
    return (s.count("G") + s.count("C")) / n if n else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--tokenizer", default="DOEJGI/GenomeOcean-4B")
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--seq-length", type=int, default=1024)
    ap.add_argument("--max-length", type=int, default=128)
    ap.add_argument("--batch-size", type=int, default=64)
    args = ap.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    rows = load_rbs_positives(args.data_dir)
    seqs = [r["seq"] for r in rows]
    print(f"RBS positives: {len(seqs):,}", flush=True)

    init_distributed(); init_parallel_state()
    print(f"loading frozen MoE from {args.checkpoint} ...", flush=True)
    model, tokenizer = load_moe_model(args)
    device = next(model.parameters()).device
    print(f"  model on {device}; extracting zonal features ...", flush=True)

    emb, rout = extract_zonal(model, tokenizer, seqs, args.batch_size, args.max_length, device)

    save = {"ids": np.array([r["id"] for r in rows]),
            "organism": np.array([r["organism"] for r in rows]),
            "phylum": np.array([r["phylum"] for r in rows]),
            "label": np.array([r["label"] for r in rows]),
            "split": np.array([r["split"] for r in rows]),
            "gc": np.array([gc_fraction(s) for s in seqs], dtype=np.float32)}
    for z in emb:
        save[f"emb_{z}"] = emb[z]
        save[f"routing_{z}"] = rout[z]
    out = os.path.join(args.output_dir, "phase0_features_rbs_zonal.npz")
    np.savez_compressed(out, **save)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
