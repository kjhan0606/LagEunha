#!/usr/bin/env bash
#SBATCH --job-name=LagQ2ENRateA10
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:A10:1
#SBATCH --mem=12G
#SBATCH --time=00:05:00
#SBATCH --nice=10000
set -euo pipefail
module purge
module load cuda/13.0.2
sha256sum -c inputs.sha256
test ! -e aligned.log
test ! -e jitter.log
nvidia-smi -L
exported=/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_w0_export_v1
./gfs_energy_neutral_full_check "$exported/aligned.bin" > aligned.log 2>&1
./gfs_energy_neutral_full_check "$exported/jitter.bin" > jitter.log 2>&1
cat aligned.log jitter.log
