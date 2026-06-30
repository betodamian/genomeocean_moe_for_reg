#!/bin/bash
#SBATCH --job-name=phase0_prom_robust
#SBATCH --account=pc_jgiga
#SBATCH --partition=es1
#SBATCH --qos=es_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:A40:1
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/global/scratch/users/bdamian6657/gomoe/logs/%x_%j.out
#SBATCH --error=/global/scratch/users/bdamian6657/gomoe/logs/%x_%j.out
#
# Phase-0 promoter robustness: (1) zonal per-position localization + (2) in-silico
# motif mutagenesis. Confounds: is the promoter signal AT the promoter motif and
# CAUSALLY dependent on it (vs a co-occurring element)? Mirrors submit_phase0_zonal.sh.
#   sbatch pipeline/submit_phase0_promoter_robustness.sh
#
set -euo pipefail

REPO=/global/scratch/users/bdamian6657/gomoe
SIF=/global/scratch/users/junhohong2028/genomeocean_moe/nemo_25.09.sif
HFC=/global/scratch/users/junhohong2028/genomeocean_moe/.cache/huggingface
CKPT=/global/scratch/users/junhohong2028/genomeocean_moe/outputs/lr_match_dropout05/stage1/step-step=49999-last
CONFIG=$REPO/configs/config_100m_moe.yaml
DATA=$REPO/data/phase0
RESULTS=$REPO/experiments/results/phase0

export APPTAINERENV_HF_HOME=$HFC
export APPTAINERENV_HF_HUB_CACHE=$HFC
export APPTAINERENV_HF_HUB_OFFLINE=1
export APPTAINERENV_TRANSFORMERS_OFFLINE=1
export APPTAINERENV_TORCH_HOME=$REPO/.cache/torch
export APPTAINERENV_PYTHONPATH=$REPO/.pylibs:/global/scratch/users/junhohong2028/genomeocean_moe/.local/lib
mkdir -p $REPO/logs $REPO/.cache/torch $RESULTS
cd $REPO

echo "=== Phase-0 promoter robustness | start $(date) | node $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

echo "--- test #1: zonal localization ---"
srun apptainer exec --nv --bind /global/scratch $SIF \
    python phase0_extract_promoter_zonal.py \
        --checkpoint "$CKPT" --config "$CONFIG" \
        --data-dir "$DATA" --output-dir "$RESULTS" \
        --batch-size 64 --max-length 128

echo "--- test #2: motif mutagenesis (val positives) ---"
srun apptainer exec --nv --bind /global/scratch $SIF \
    python phase0_extract_promoter_mutagenesis.py \
        --checkpoint "$CKPT" --config "$CONFIG" \
        --data-dir "$DATA" --output-dir "$RESULTS" \
        --split val --batch-size 64 --max-length 128

echo "=== done $(date) ==="
ls -la $RESULTS/phase0_features_promoters_zonal.npz $RESULTS/phase0_features_promoters_mutagenesis.npz
