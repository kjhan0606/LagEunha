#!/usr/bin/env bash
#SBATCH --job-name=LagQ2ENGuardBuild
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --nice=10000
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
sha256sum -c inputs.sha256
run_dir=$PWD
source_tar=/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_native_build_v1/source.tar
build_dir=${SLURM_TMPDIR:-/tmp}/lag_q2_energy_neutral_guard_${SLURM_JOB_ID}
mkdir -p "$build_dir"
tar -xf "$source_tar" -C "$build_dir"
/home/kjhan/miniconda3/bin/python3 apply_kepler_gpu_energy_neutral_guard_20261027.py \
    "$build_dir" --receipt "$run_dir/guard_receipt.json"
cd "$build_dir"
TIMEFORMAT='[BUILD_TIME] wall=%R'
time make -j2 USE_CUDA=1 GFS_VARPRESS_W0=1 > "$run_dir/build.log" 2>&1
test -s eunha2.exe
cp eunha2.exe "$run_dir/eunha2_energy_neutral_guard.exe"
sha256sum "$run_dir/eunha2_energy_neutral_guard.exe" > "$run_dir/binary.sha256"
tail -n 12 "$run_dir/build.log"
