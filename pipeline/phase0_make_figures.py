#!/usr/bin/env python3
"""
Make Phase-0 presentation figures (3 PNGs) from the committed result JSONs.
Run: .venv/bin/python pipeline/phase0_make_figures.py
Output: experiments/phase0/figures/*.png
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PH   = os.path.join(ROOT, "experiments/phase0")
FIG  = os.path.join(PH, "figures")
os.makedirs(FIG, exist_ok=True)

smoke = json.load(open(os.path.join(PH, "phase0_results.json")))
zonal = json.load(open(os.path.join(PH, "phase0_sd_zonal_results.json")))

# pretty task names (order: 4 passing then the SD one)
TASKS = [
    ("promoter_vs_intergenic",     "Promoter\nvs intergenic"),
    ("rbs_TIS_vs_intergenic",      "Gene start (RBS-TIS)\nvs intergenic"),
    ("rho_t1_vs_intergenic",       "Rho terminator\nvs intergenic"),
    ("rho_t1_vs_intrinsic_ecoli",  "Rho vs intrinsic\nterminator (E. coli)"),
    ("rbs_SD_vs_UNSD",             "SD motif\n(SD vs UNSD)"),
]
BLUE, GREY, GREEN, ORANGE = "#2c6fbb", "#b0b0b0", "#3a9d6e", "#e08a1e"

# ── Figure 1: go/no-go — model (routing_concat) vs GC baseline ────────────────
fig, ax = plt.subplots(figsize=(9, 5))
labels = [t[1] for t in TASKS]
x = np.arange(len(TASKS)); w = 0.38
rc  = [smoke[t]["views"]["routing_concat"]["mcc"] for t, _ in TASKS]
rc_lo = [smoke[t]["views"]["routing_concat"]["mcc"] - smoke[t]["views"]["routing_concat"]["mcc_ci"][0] for t, _ in TASKS]
rc_hi = [smoke[t]["views"]["routing_concat"]["mcc_ci"][1] - smoke[t]["views"]["routing_concat"]["mcc"] for t, _ in TASKS]
gc  = [smoke[t]["views"]["gc_only"]["mcc"] for t, _ in TASKS]
b1 = ax.bar(x - w/2, rc, w, yerr=[rc_lo, rc_hi], capsize=3, color=BLUE, label="GenomeOcean-MoE (routing+embedding)")
b2 = ax.bar(x + w/2, gc, w, color=GREY, label="GC-content baseline (chance)")
for i, (t, _) in enumerate(TASKS):
    ok = smoke[t]["gate_pass"]
    ax.annotate("PASS" if ok else "FAIL", (i - w/2, rc[i] + rc_hi[i] + 0.03),
                ha="center", fontsize=9, fontweight="bold",
                color=GREEN if ok else ORANGE)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8.5)
ax.set_ylabel("Discrimination skill (MCC)\n0 = chance, 1 = perfect")
ax.set_ylim(-0.1, 1.2)
ax.axhline(0, color="k", lw=0.6)
ax.set_title("Phase-0 go/no-go: the frozen model already separates 4 of 5 regulatory\n"
             "tasks from same-context decoys, far above a GC baseline", fontsize=11)
ax.legend(loc="upper right", fontsize=8.5, framealpha=0.95)
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig1_go_nogo.png"), dpi=150); plt.close()

# ── Figure 2: routing vs embedding (MoE-necessity, P1) ────────────────────────
fig, ax = plt.subplots(figsize=(9, 5))
w = 0.27
emb = [smoke[t]["views"]["embedding_only"]["mcc"] for t, _ in TASKS]
ro  = [smoke[t]["views"]["routing_only"]["mcc"] for t, _ in TASKS]
rc  = [smoke[t]["views"]["routing_concat"]["mcc"] for t, _ in TASKS]
ax.bar(x - w, emb, w, color=GREY,  label="Embedding only (a dense model also has this)")
ax.bar(x,     ro,  w, color=ORANGE, label="Routing only (MoE-specific channel)")
ax.bar(x + w, rc,  w, color=BLUE,  label="Routing + embedding (full)")
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8.5)
ax.set_ylabel("Discrimination skill (MCC)")
ax.set_ylim(-0.05, 1.0); ax.axhline(0, color="k", lw=0.6)
ax.set_title("The MoE routing channel carries real signal beyond the embedding\n"
             "(routing > embedding in the 4 passing tasks — the MoE-necessity prediction)",
             fontsize=11)
ax.legend(loc="upper right", fontsize=8.5)
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig2_routing_vs_embedding.png"), dpi=150); plt.close()

# ── Figure 3: the SD test — signal is there, but the router is blind ───────────
m = zonal["macro"]
order = [("kmer4_sd", "Raw DNA k-mers\n(SD region)", GREEN),
         ("emb_sd", "Model embedding\n(at SD, zoomed in)", GREY),
         ("emb_full", "Model embedding\n(whole window)", GREY),
         ("routing_full", "Routing\n(whole window)", ORANGE),
         ("routing_sd", "Routing\n(at SD, zoomed in)", ORANGE)]
fig, ax = plt.subplots(figsize=(9, 5))
vals = [m[k] for k, _, _ in order]
cols = [c for _, _, c in order]
bars = ax.bar(range(len(order)), vals, color=cols)
for i, v in enumerate(vals):
    ax.annotate(f"{v:.2f}", (i, v + 0.02), ha="center", fontsize=10, fontweight="bold")
ax.set_xticks(range(len(order))); ax.set_xticklabels([l for _, l, _ in order], fontsize=8.5)
ax.set_ylabel("SD signal captured\n(Spearman ρ vs measured SD binding energy)")
ax.set_ylim(0, 1.0)
ax.set_title("Why the SD sub-task failed: the signal is clearly in the DNA (k-mers, 0.87),\n"
             "but the MoE routing channel is blind to it (0.01) — even zoomed onto the SD",
             fontsize=11)
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig3_sd_router_blind.png"), dpi=150); plt.close()

print("wrote:")
for f in ["fig1_go_nogo.png", "fig2_routing_vs_embedding.png", "fig3_sd_router_blind.png"]:
    print(" ", os.path.join("experiments/phase0/figures", f))
