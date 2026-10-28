#!/usr/bin/env bash
#SBATCH --job-name=LagQ2WdotSmoke
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --nice=10000
set -euo pipefail
module purge
module load cuda/13.0.2
sha256sum -c inputs.sha256
nvidia-smi -L
./gfs_volume_rate_probe_smoke > result.log 2>&1
cat result.log
