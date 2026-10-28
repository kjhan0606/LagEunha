#!/usr/bin/env bash
#SBATCH --job-name=LagQ2ZeroPBuild
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
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
sha256sum -c inputs.sha256
run_dir=$PWD
build_dir=${SLURM_TMPDIR:-/tmp}/lag_q2_zero_pressure_${SLURM_JOB_ID}
mkdir -p "$build_dir"
tar -xf source.tar -C "$build_dir"
python3 patch_guard.py "$build_dir" > patch_guard_receipt.json
python3 patch_iterations.py "$build_dir" > patch_iterations_receipt.json
python3 patch_zero.py "$build_dir" > patch_zero_receipt.json
g++ -std=c++14 -O2 -Wall -Wextra -Werror -I "$build_dir/Exam" \
    test_projection.cpp -o "$build_dir/test_projection"
"$build_dir/test_projection" > test_projection.log
cd "$build_dir"
TIMEFORMAT='[BUILD_TIME] wall=%R'
time make -j2 USE_CUDA=1 GFS_VARPRESS_W0=1 > "$run_dir/build.log" 2>&1
test -s eunha2.exe
cp eunha2.exe "$run_dir/eunha2_zero_pressure_v4.exe"
sha256sum "$run_dir/eunha2_zero_pressure_v4.exe" > "$run_dir/binary.sha256"
tail -n 12 "$run_dir/build.log"
