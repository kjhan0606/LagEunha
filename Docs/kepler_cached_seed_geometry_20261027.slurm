#!/usr/bin/env bash
#SBATCH --job-name=LagQ2CachedSeedGeom
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
root=/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_seedcert_w0_weighted_pair_v1
/home/kjhan/miniconda3/bin/python3 audit_kepler_cached_seed_geometry_20261027.py "$root/w0/baseline" --output w0.json
/home/kjhan/miniconda3/bin/python3 audit_kepler_cached_seed_geometry_20261027.py "$root/weighted/baseline" --output weighted.json
