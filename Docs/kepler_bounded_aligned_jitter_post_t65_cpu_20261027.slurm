#!/usr/bin/env bash
#SBATCH --job-name=LagQ2Cmp65
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
work=$PWD
root=/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027
aligned="$root/to65_near32_a10_v1/audit_t65.json"
jitter="$root/jitter005_to65_a10_v1/audit_t65.json"
test -s "$aligned"
test -s "$jitter"
sha256sum -c prestaged.sha256
/home/kjhan/miniconda3/bin/python3 "$work/audit_kepler_bounded_radial_tails_20261027.py" \
    "$aligned" --output aligned_radial_tails.json > aligned_radial_tails.log
/home/kjhan/miniconda3/bin/python3 "$work/audit_kepler_bounded_radial_tails_20261027.py" \
    "$jitter" --output jitter_radial_tails.json > jitter_radial_tails.log
/home/kjhan/miniconda3/bin/python3 "$work/audit_kepler_bounded_aligned_jitter_pair_20261027.py" \
    "$aligned" "$jitter" aligned_radial_tails.json jitter_radial_tails.json \
    --output pair_t65.json > pair_t65.log
/home/kjhan/miniconda3/bin/python3 "$work/plot_kepler_bounded_aligned_jitter_density_20261027.py" \
    "$aligned" "$jitter" --output density_t65.png > density_t65.log
