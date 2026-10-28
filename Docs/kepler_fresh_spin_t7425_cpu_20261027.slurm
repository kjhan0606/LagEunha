#!/usr/bin/env bash
#SBATCH --job-name=LagQ2Spin7425
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:0
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --nice=200
set -euo pipefail
module purge
module load gnu13/13.2.0
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
sha256sum -c inputs.sha256
python=/home/kjhan/miniconda3/bin/python3
"$python" test_kepler_affine_spin_energy_20261027.py
"$python" audit_kepler_affine_spin_energy_20261027.py weighted_receipt.json --thermal-risk --output weighted.json
