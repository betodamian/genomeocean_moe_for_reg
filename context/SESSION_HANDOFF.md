# Session handoff — 2026-06-29 (Week 1 complete + Phase-0 GO)

Pick-up doc for continuing this project in a fresh session. Read `context/research_plan.md`
(the v2.8 plan) first; this file says what's **done**, the **state**, and **what's next**.

---

## TL;DR

- **Week-1 data pipeline is complete** for all 3 elements (promoters, RBS, Rho): windows
  built, MMseqs2 ≤60% clustering done, census gates PASS, 80/20 splits + LOO folds committed.
- **Phase-0 ceiling benchmark = GO.** 4/5 same-context tasks clear the GC-matched gate;
  routing beats embedding in the 4 passing tasks (the MoE-necessity direction, P1).
- **The 1 failure (RBS-SD) was fully diagnosed and resolved** as an honest negative:
  the MoE router does not specialize for the short Shine-Dalgarno motif (non-specialization,
  not pooling dilution). RBS-TIS detection itself is strong and is the headline RBS result.
- Everything is committed + pushed to GitHub `betodamian/genomeocean_moe_for_reg` (main),
  authored as the user (no Claude co-author — keep it that way).

---

## What was done this session (in order)

1. **S. aureus Ribo-RET integration finished.** Re-ran `pipeline/parse_label_tables.py`
   (fixed `SAUR_EXSD` organism label → `saureus_HG001`). `peakcall_saureus_tis.py`
   (written previously) gives 2,122 TIS. Updated `data/MANIFEST.md` (HG001 = 9th genome).

2. **RBS windows + split.** `pipeline/build_rbs_windows.py` → 20,696 TIS windows (300 bp,
   TIS centered at index 150) across 7 organisms, labeled positive_SD / positive_UNSD by an
   upstream motif scan, + 20,696 Tier-1 intergenic decoys. `pipeline/cluster_rbs.py` →
   17,695 ≤60% clusters, census **PASS** (4 bacterial phyla + archaea). Artifacts in
   `splits/rbs/` (clusters.tsv, split_80_20.tsv, folds_loo.tsv, census.md). NOTE: E. coli
   K-12 + BL21 are co-held-out in one LOO fold (near-identical core genes → leakage).

3. **Rho windows + split.** `pipeline/peakcall_ecoli_rho_termseq.py` → 151 bp-resolution
   E. coli Rho sites from Term-seq (GSE109766) ∩ BCM-BST concordance (avoids centering on
   huge BST regions). `pipeline/build_rho_windows.py` → 5,682 positives (T1 in-vivo:
   151 E. coli + 742 MTB; T2 in-vitro: 4,789) + TERMITe intrinsic decoys + intergenic.
   `pipeline/cluster_rho.py` → census **PASS** (2 in-vivo phyla, honest 2-organism scope).
   Artifacts in `splits/rho/`.

4. **Phase-0 smoke test (the gate).** Two-stage:
   - `pipeline/phase0_extract_features.py` (runs on LBL cluster, NeMo container) — one frozen
     forward pass per window on GenomeOcean-MoE pilot4; saves embedding(768)+routing(96)+
     gc+kmer4 per window → `phase0_features_<element>.npz`. **Must call init_distributed()
     + init_parallel_state() before load_moe_model** (Megatron parallel state).
   - `pipeline/phase0_smoke_test.py` (CPU) — 5 pre-registered same-context tasks, linear
     probes across 5 feature views, bootstrap 95% CIs, gate = routing_concat MCC lower-CI >
     gc_only MCC upper-CI.
   - `pipeline/submit_phase0.sh` — A40 SLURM job.
   - **Result: GO.** promoter 0.63, RBS-TIS 0.81, Rho-vs-intergenic 0.70, Rho-vs-intrinsic
     0.78 all PASS; RBS-SD-vs-UNSD FAIL.

5. **RBS-SD investigation (data → method).**
   - Found the SD labels were **circular** (regex motif-presence) and **cross-organism
     GC-confounded**. Fix: `pipeline/compute_sd_deltaG.py` relabels SD strength as the
     hybridization ΔG between each gene's upstream and that organism's own 16S 3' anti-SD
     (ViennaRNA), per organism → `data/datasets/rbs/sd_deltaG.tsv`.
   - `pipeline/phase0_sd_reanalysis.py` (within-organism ΔG regression): confound removed
     (gc ρ 0.28→0.02) but model still weak, routing < embedding.
   - **Decisive per-position test:** `pipeline/phase0_extract_sd_zonal.py` (cluster) captures
     per-token routing pooled over only the SD-zone tokens (token→bp coords from BPE strings,
     0/20,696 mismatches). `pipeline/phase0_sd_zonal_analysis.py` →
     **verdict (b) non-specialization**: routing_sd 0.01 ≈ routing_full 0.06 (no dilution),
     routing_sd 0.01 ≪ emb_sd 0.40 (P1 fails at SD zone), kmer4_sd 0.87 (signal IS in DNA).
     `pipeline/submit_phase0_zonal.sh` is the job.

6. **Figures.** `pipeline/phase0_make_figures.py` → `experiments/phase0/figures/`:
   fig1_go_nogo, fig2_routing_vs_embedding, fig3_sd_router_blind.

Reports/artifacts: `experiments/phase0/PHASE0_FINDINGS.md` (full writeup), `phase0_report.md`,
`phase0_results.json`, `phase0_sd_report.md`, `phase0_sd_zonal_report.md` (+ JSONs).

---

## Cluster operational notes (LBL)

- Access: `ssh lrc` (host alias; user `bdamian6657`). **scp fails** (mux/login) — transfer
  with `tar czf - ... | ssh lrc 'cd DIR && tar xzf -'` and pull with the reverse.
- Working dir on cluster: `/global/scratch/users/bdamian6657/gomoe`
  (has `eval/` reusable scripts, `configs/config_100m_moe.yaml`, `pipeline/`, `data/phase0/`).
- Checkpoint (frozen MoE, pilot4): `/global/scratch/users/junhohong2028/genomeocean_moe/outputs/lr_match_dropout05/stage1/step-step=49999-last` (read-only, accessible).
- Container: `/global/scratch/users/junhohong2028/genomeocean_moe/nemo_25.09.sif`.
- Reusable eval functions live in `gomoe/eval/`: `load_moe_model`, `init_distributed`,
  `init_parallel_state`, `build_mixtral_config` (eval_utils.py); `extract_features_moe`,
  `extract_embeddings_moe`, `run_probe` (eval_downstream.py); `build_token_coordinate_map`
  (eval_expert_specialization.py).
- SLURM: `--account=pc_jgiga --partition=es1 --qos=es_normal --gres=gpu:A40:1`. Jobs take
  ~6 min for full extraction. Always `set -euo pipefail` in submit scripts (so a Python
  crash surfaces as FAILED, not a false COMPLETED).
- Features: 864-d = 768 embedding + 96 routing (12 layers × 8 experts). Tokenizer
  `DOEJGI/GenomeOcean-4B`, ~4.89 bp/token.

## Local environment

- `.venv` has: numpy, scipy, scikit-learn, matplotlib, ViennaRNA, openpyxl, xlrd.
- NPZ feature files (`data/datasets/phase0/*.npz`, ~270 MB) are **gitignored** — regenerable
  via the submit scripts. `data/datasets/` and `*.npz` are gitignored; `splits/` IS committed.

## Git workflow (REQUIRED)

Commit authored as the user only — **never** add a Claude co-author trailer. Push every
significant change to `origin/main`. (See memory `git-workflow.md`.)

---

## Phase-0 results (numbers)

| task | routing_concat MCC | gc baseline | routing_only | embedding_only | gate |
|---|---|---|---|---|---|
| promoter vs intergenic | 0.63 | −0.02 | 0.54 | 0.44 | PASS |
| RBS-TIS vs intergenic | 0.81 | 0.02 | 0.71 | 0.61 | PASS |
| Rho vs intergenic | 0.70 | 0.55 | 0.67 | 0.48 | PASS |
| Rho vs intrinsic (E. coli) | 0.78 | 0.10 | 0.82 | 0.41 | PASS |
| RBS SD vs UNSD | 0.34 | 0.28 | 0.27 | 0.34 | FAIL → resolved (non-specialization) |

SD zonal (within-organism ΔG regression, Spearman ρ macro): kmer4_sd 0.87, emb_sd 0.40,
emb_full 0.15, routing_full 0.06, routing_sd 0.01.

---

## Next steps (research_plan §16)

1. **Week 3 — detection vs baselines (P2).** Run classical tools (BPROM, Prodigal-RBS,
   RhoTermPredict) and dense gLMs (NT, DNABERT-2, Evo 2, ProkBERT, GO-dense) on the SAME
   held-out windows; compare AUPRC + boundary-F1. Lead with the 4 passing tasks.
2. **Week 4 — MoE necessity (P3/P4) + generalization (P8).** Expert-ablation DiD on the
   passing elements (extend `eval/eval_expert_ablation.py`); leave-one-organism-out
   cross-genome folds (already in `splits/*/folds_loo.tsv`) for the unseen-split test.
3. **Optional enrichment (deferred):** E. coli ΔaSD Ribo-seq (GSE135906, WIG pulled),
   MTB dRNA-seq leaderless (SRP028740) → more UNSD/leaderless labels.
4. **Optional (P10):** within-token conditional-MI control on the passing elements
   (rule out "trivial" motif-token specialization).

Headline framing for write-up: GO; routing carries regulatory signal beyond the embedding
for promoters / gene-starts / Rho; RBS-TIS is the strong RBS result; RBS-SD is a clean,
interpretable negative (router specializes for structural context, not a 6-bp motif).
