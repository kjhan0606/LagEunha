#!/usr/bin/env bash
#SBATCH --job-name=LagQ2Slip65
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --nice=10000
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
sha256sum -c source.sha256
variant="${1:?Expected aligned or jitter}"
case "$variant" in
    aligned|jitter) ;;
    *) echo "Invalid variant: $variant" >&2; exit 2 ;;
esac
/home/kjhan/miniconda3/bin/python3 audit_kepler_bounded_weighted_face_slip_t65_20261027.py \
    "${variant}_receipt.json" "${variant}_weight.npz" "${variant}_rate.npz" \
    --output "${variant}_slip.json" > "${variant}_slip.log"
