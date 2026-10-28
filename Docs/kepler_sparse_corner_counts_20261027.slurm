#!/usr/bin/env bash
#SBATCH --job-name=KeplerCornerCount
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:20:00
set -euo pipefail
module purge
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/home/kjhan/miniconda3/bin/python3 \
  /home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs/audit_kepler_sparse_corner_counts_20261027.py \
  "$@"
