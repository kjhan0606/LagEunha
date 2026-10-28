#!/usr/bin/env bash
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
sha256sum -c inputs.sha256
run_dir=$PWD
build_dir=/tmp/lag_q2_limiter_only_v7_manual_001
test ! -e "$build_dir"
mkdir "$build_dir"
tar -xf source.tar -C "$build_dir"
python3 patch_guard.py "$build_dir" > patch_guard_receipt.json
python3 patch_iterations.py "$build_dir" > patch_iterations_receipt.json
python3 patch_limiter.py "$build_dir" > patch_limiter_receipt.json
cd "$build_dir"
TIMEFORMAT='[BUILD_TIME] wall=%R'
time make -j2 USE_CUDA=1 GFS_VARPRESS_W0=1 > "$run_dir/build.log" 2>&1
test -s eunha2.exe
cp eunha2.exe "$run_dir/eunha2_energy_neutral_limiter_v7.exe"
sha256sum "$run_dir/eunha2_energy_neutral_limiter_v7.exe" > "$run_dir/binary.sha256"
tail -n 14 "$run_dir/build.log"
