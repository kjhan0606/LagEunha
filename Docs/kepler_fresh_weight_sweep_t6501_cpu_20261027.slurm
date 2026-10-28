#!/usr/bin/env bash
#SBATCH --job-name=LagQ2WeightSweep65
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:0
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --nice=200
set -euo pipefail
module purge
module load gnu13/13.2.0
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
sha256sum -c inputs.sha256
/home/kjhan/miniconda3/bin/python3 audit_kepler_fresh_weight_sweep_20261027.py /gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_capacity64_v5_fresh_localwork_t6501_ada_v1/audit_t6501.json --output weight_sweep.json
