#!/usr/bin/env bash
#SBATCH --job-name=LagPolyMomentCPU
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:0
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --nice=10000
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
sha256sum -c inputs.sha256
g++ -std=c++14 -O2 -Wall -Wextra -Werror -shared -fPIC -I . \
    gfs_polygon_moment_cpu_bridge_20261027.cpp -o libpolygon_moment.so
/home/kjhan/miniconda3/bin/python3 audit_kepler_polygon_moment_native_faces_20261027.py \
    w0_jitter_t65.json --library ./libpolygon_moment.so --output w0_jitter_t65_moments.json > w0.log
/home/kjhan/miniconda3/bin/python3 audit_kepler_polygon_moment_native_faces_20261027.py \
    weighted_jitter_t10.json --expected-failed-check disk_density_equilibrium \
    --library ./libpolygon_moment.so --output weighted_jitter_t10_moments.json > weighted.log
