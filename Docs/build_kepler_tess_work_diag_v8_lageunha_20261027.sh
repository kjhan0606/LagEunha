#!/usr/bin/env bash
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
sha256sum -c inputs.sha256
run_dir=$PWD
for mode in seed5 seed3; do
    build_dir=$(mktemp -d "/tmp/lag_q2_tess_work_diag_v8_${mode}_XXXXXXXX")
    tar -xf source.tar -C "$build_dir"
    python3 patch_guard.py "$build_dir" > "$run_dir/patch_guard_${mode}.json"
    python3 patch_iterations.py "$build_dir" > "$run_dir/patch_iterations_${mode}.json"
    if [[ "$mode" == seed3 ]]; then
        python3 patch_diag.py "$build_dir" --seed3 > "$run_dir/patch_${mode}.json"
    else
        python3 patch_diag.py "$build_dir" > "$run_dir/patch_${mode}.json"
    fi
    cd "$build_dir"
    TIMEFORMAT='[BUILD_TIME] wall=%R'
    time make -j2 USE_CUDA=1 GFS_VARPRESS_W0=1 > "$run_dir/build_${mode}.log" 2>&1
    test -s eunha2.exe
    cp eunha2.exe "$run_dir/eunha2_tess_work_${mode}.exe"
    sha256sum "$run_dir/eunha2_tess_work_${mode}.exe" > "$run_dir/binary_${mode}.sha256"
    cd "$run_dir"
    tail -n 5 "build_${mode}.log"
done
