#!/usr/bin/env bash
#SBATCH --job-name=LagQ2ForceSeries
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --nice=10000
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
sha256sum -c prestaged.sha256
root=/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027
for case in aligned jitter; do
    if [ "$case" = aligned ]; then
        run="$root/to10_near32_a10_v1"
        extra=()
    else
        run="$root/jitter005_to10_a10_v1"
        extra=(--expected-failed-check disk_ballistic_residual)
    fi
    for idx in 000053 000478 000904 001329 001754 001806; do
        /home/kjhan/miniconda3/bin/python3 audit_kepler_bounded_pressure_force_t10_20261027.py \
            "$run/audit_t10.json" --checkpoint "$run/khout.$idx.dat" \
            "${extra[@]}" --output "${case}_${idx}.json" > "${case}_${idx}.log"
    done
done
/home/kjhan/miniconda3/bin/python3 summarize_kepler_bounded_pressure_force_series_20261027.py \
    > pressure_force_series_t10.log
