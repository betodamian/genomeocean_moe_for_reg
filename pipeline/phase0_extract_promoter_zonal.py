#!/usr/bin/env python3
"""
Phase-0 promoter ROBUSTNESS test #1 — per-position (zonal) localization
(research_plan §4 / §9b confound control). Runs ON the LBL cluster, frozen MoE.

Question: the Phase-0 promoter task (promoter vs same-context intergenic, MCC 0.63)
shows the model separates promoter windows — but is it keying on the PROMOTER MOTIF,
or on some co-occurring "mystery" feature elsewhere in the 300-bp window? If the
discriminative signal is carried by tokens AT the core-promoter region (−35/−10/TSS)
and NOT by downstream / far-upstream tokens, the model is using the promoter, not a
co-located confound.

Like phase0_extract_sd_zonal.py, this captures PER-TOKEN routing + hidden states and
pools them within ZONES (token→bp coords reconstructed from BPE strings), but for
BOTH promoter positives AND their intergenic decoys (the classification needs both).

Window geometry (build_promoter_windows.py): 300 bp, promoter strand 5'→3',
TSS at window index 150. Zones (bp index in the window):
  minus35 = [112,124)  (−38..−26)   bacterial −35 hexamer
  minus10 = [136,146)  (−14..−4)    bacterial −10 hexamer
  core    = [110,150)  (−40..−1)    core promoter (covers bacterial −35/−10 AND
                                     archaeal TATA/BRE — organism-agnostic)
  tss     = [146,156)  (−4..+6)     transcription start
  down    = [156,300)  (+6..+150)   transcript body (NO promoter motif → control)
  farup   = [0,110)    (−150..−40)  far upstream (mostly non-motif → control)
  full    = all content tokens (reproduces the original mean-pool; alignment control)

Decision (in phase0_promoter_zonal_analysis.py):
  MCC(core) ≫ MCC(down) and MCC(core) ≫ MCC(farup)  → signal localized to the
    promoter region (NOT a mystery element elsewhere).
  routing_core > emb_core                            → P1 (MoE channel) at the motif.

Output: <out>/phase0_features_promoters_zonal.npz
Run inside container — see submit_phase0_promoter_robustness.sh.
"""
import argparse, csv, os, sys
import numpy as np
import torch

EVAL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "eval")
sys.path.insert(0, EVAL_DIR)
from eval_utils import init_distributed, init_parallel_state, load_moe_model  # noqa: E402

ZONES = {"minus35": (112, 124), "minus10": (136, 146), "core": (110, 150),
         "tss": (146, 156), "down": (156, 300), "farup": (0, 110)}  # full handled separately


def load_promoter_windows(data_dir):
    all_tsv   = os.path.join(data_dir, "promoters_ALL.tsv")
    split_tsv = os.path.join(data_dir, "promoters_split_80_20.tsv")
    split_of = {}
    with open(split_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            split_of[r["window_id"]] = r["split"]
    rows = []
    with open(all_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(dict(id=r["id"], organism=r["organism"], phylum=r["phylum"],
                             label=r["label"], split=split_of.get(r["id"], "none"),
                             seq=r["window_seq"]))
    return rows


def token_bp_coords(token_strings):
    coords, bp = [], 0
    for ts in token_strings:
        clean = ts.replace("▁", "").replace("Ġ", "").strip().upper()
        clean = "".join(c for c in clean if c in "ACGTN")
        coords.append((bp, bp + len(clean))); bp += len(clean)
    return coords


def zone_mask(coords, z0, z1):
    return np.array([(s < z1 and e > z0) for (s, e) in coords], dtype=bool)


def extract_zonal(model, tokenizer, sequences, batch_size, max_length, device):
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
            hidden = captured["hidden"].transpose(0, 1)
            parts = []
            for li in range(num_layers):
                parts.append(router_captured[li].reshape(B, S, num_experts) if li in router_captured
                             else torch.zeros(B, S, num_experts, device=device))
            routing = torch.cat(parts, dim=-1).cpu().float().numpy()
            hidden_np = hidden.cpu().float().numpy()
            am = enc["attention_mask"].numpy().astype(bool)
            ids_np = input_ids.cpu().numpy()
            for b in range(B):
                content = am[b].copy()
                for sid in special:
                    content &= (ids_np[b] != sid)
                cidx = np.where(content)[0]
                toks = tokenizer.convert_ids_to_tokens(ids_np[b][cidx].tolist())
                coords = token_bp_coords(toks)
                emb_t = hidden_np[b][cidx]
                rout_t = routing[b][cidx]
                if len(coords) != emb_t.shape[0]:
                    n_mismatch += 1
                out_emb["full"].append(emb_t.mean(0) if len(emb_t) else np.zeros(768, np.float32))
                out_rout["full"].append(rout_t.mean(0) if len(rout_t) else np.zeros(96, np.float32))
                for z, (z0, z1) in ZONES.items():
                    m = zone_mask(coords, z0, z1)
                    if m.any():
                        out_emb[z].append(emb_t[m].mean(0)); out_rout[z].append(rout_t[m].mean(0))
                    else:
                        out_emb[z].append(np.zeros(768, np.float32)); out_rout[z].append(np.zeros(96, np.float32))
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

    rows = load_promoter_windows(args.data_dir)
    seqs = [r["seq"] for r in rows]
    print(f"promoter windows (pos+decoy): {len(seqs):,}", flush=True)

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
        save[f"emb_{z}"] = emb[z]; save[f"routing_{z}"] = rout[z]
    out = os.path.join(args.output_dir, "phase0_features_promoters_zonal.npz")
    np.savez_compressed(out, **save)
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
