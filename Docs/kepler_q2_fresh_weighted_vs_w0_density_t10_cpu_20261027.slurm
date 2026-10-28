#!/usr/bin/env bash
#SBATCH --job-name=LagQ2FreshDensity10
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:0
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --nice=10000
set -euo pipefail
module purge
module load gnu13/13.2.0
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/mpl_${SLURM_JOB_ID}"
sha256sum -c inputs.sha256
root=/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027
python=/home/kjhan/miniconda3/bin/python3
"$python" plot_kepler_q2_weighted_vs_w0_density_20261027.py \
  "$root/jitter005_to10_a10_v1/audit_t10.json" \
  "$root/gpu_energy_neutral_capacity64_v5_fresh_to10_a10_v1/audit_t10.json" \
  --w0-expected-failed-check disk_ballistic_residual \
  --output density.png --pixels 768 --relative-max .005
sha256sum density.png > output.sha256
