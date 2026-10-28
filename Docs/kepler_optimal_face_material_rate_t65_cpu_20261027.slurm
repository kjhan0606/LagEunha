#!/usr/bin/env bash
#SBATCH --job-name=LagQ2FaceOpt
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:0
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --nice=10000
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
sha256sum -c inputs.sha256
for case in aligned jitter; do
  /home/kjhan/miniconda3/bin/python3 audit_kepler_optimal_face_material_rate_20261027.py \
    "$case" "${case}_receipt.json" --output "${case}_optimal_face_rate.json" > "${case}.log"
done
