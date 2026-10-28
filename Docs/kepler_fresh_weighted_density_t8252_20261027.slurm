#!/usr/bin/env bash
#SBATCH --job-name=LagQ2Density8252
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:0
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --nice=200
set -euo pipefail
module purge
module load gnu13/13.2.0
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/mpl_${SLURM_JOB_ID}"
sha256sum -c inputs.sha256
compare=/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/fresh_compare_step904_t8252_cpu_v1
/home/kjhan/miniconda3/bin/python3 plot_kepler_q2_weighted_vs_w0_density_20261027.py "$compare/w0_receipt.json" "$compare/weighted_receipt.json" --max-time-difference 1e-6 --output density.png --pixels 768 --relative-max .005
sha256sum density.png > output.sha256
