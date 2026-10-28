#!/usr/bin/env bash
#SBATCH --job-name=KeplerSparseSeries
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
module purge
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [ "$#" -eq 0 ]; then
  echo 'Expected audit_kepler_sparse_a10_series_20261027.py arguments' >&2
  exit 2
fi
/home/kjhan/miniconda3/bin/python3 \
  /home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs/audit_kepler_sparse_a10_series_20261027.py \
  "$@"
