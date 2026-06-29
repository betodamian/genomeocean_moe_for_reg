#!/bin/bash
#SBATCH --job-name=phase0_extract
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
# Phase-0 feature extraction for promoters + RBS + Rho windows on the frozen MoE.
# Mirrors submit_r7b.sh conventions (NeMo container, junho's SIF/ckpt/HF cache).
#
#   sbatch pipeline/submit_phase0.sh
#
set -uo pipefail

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

echo "=== Phase-0 extract | start $(date) | node $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

srun apptainer exec --nv --bind /global/scratch $SIF \
    python phase0_extract_features.py --element all \
        --checkpoint "$CKPT" --config "$CONFIG" \
        --data-dir "$DATA" --output-dir "$RESULTS" \
        --batch-size 64 --max-length 128

echo "=== done $(date) ==="
ls -la $RESULTS
