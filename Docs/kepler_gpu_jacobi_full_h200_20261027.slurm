#!/usr/bin/env bash
#SBATCH --job-name=LagQ2WdotFull
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --nice=10000
set -euo pipefail
module purge
module load cuda/13.0.2
sha256sum -c inputs.sha256
nvidia-smi -L
./gfs_volume_rate_probe_full_check aligned.bin > aligned.log 2>&1
./gfs_volume_rate_probe_full_check jitter.bin > jitter.log 2>&1
cat aligned.log jitter.log
