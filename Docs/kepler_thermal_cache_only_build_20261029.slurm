#!/bin/bash
# Historical attempt: grammar CPU cluster lacked CUDA 13; see THERMAL_ONLY_BUILD_ATTEMPT.json.
#SBATCH --job-name=KepThermalCacheBuild
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --chdir=/gpfs/kjhan/LagForce/Kepler_ThermalCacheOnly_20261029
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_ThermalCacheOnly_20261029/slurm-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_ThermalCacheOnly_20261029/slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=2
sha256sum -c build_source.sha256
/home/kjhan/miniconda3/bin/python3 ./build.py /gpfs/kjhan/LagForce/Kepler_ThermalCacheOnly_20261029
