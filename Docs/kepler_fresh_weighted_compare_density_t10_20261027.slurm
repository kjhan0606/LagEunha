#!/usr/bin/env bash
#SBATCH --job-name=LagQ2FreshCompare10
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:0
#SBATCH --mem=12G
#SBATCH --time=00:20:00
#SBATCH --nice=200
set -euo pipefail
module purge
module load gnu13/13.2.0
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/mpl_${SLURM_JOB_ID}"
sha256sum -c inputs.sha256
root=/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027
/home/kjhan/miniconda3/bin/python3 audit_kepler_bounded_entropy_sign_20261027.py w0_receipt.json --allow-failed-check disk_ballistic_residual --output w0_entropy.json
/home/kjhan/miniconda3/bin/python3 audit_kepler_weighted_entropy_sign_20261027.py weighted_receipt.json --output weighted_entropy.json
/home/kjhan/miniconda3/bin/python3 audit_kepler_weighted_vs_w0_fresh_20261027.py w0_receipt.json weighted_receipt.json w0_entropy.json weighted_entropy.json --allow-w0-failed-check disk_ballistic_residual --output comparison.json
/home/kjhan/miniconda3/bin/python3 plot_kepler_q2_weighted_vs_w0_density_20261027.py w0_receipt.json weighted_receipt.json --w0-expected-failed-check disk_ballistic_residual --max-time-difference 1e-6 --output density.png --pixels 768 --relative-max .005
sha256sum comparison.json w0_entropy.json weighted_entropy.json density.png > outputs.sha256
