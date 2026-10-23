#!/usr/bin/env bash
set -euo pipefail
dest=$(cd -- "$(dirname -- "$0")" && pwd)
test ! -e "$dest/continuation"
source "$dest/run_env.sh"
export I_MPI_CXX=g++
sha256sum -c "$dest/source.sha256" > "$dest/prebuild_source_integrity.txt"
cd "$dest/src"
nvcc -O2 -std=c++17 --fmad=false -DXYZDBL -DGOTPM -DUSE_CUDA \
    -gencode arch=compute_89,code=sm_89 -gencode arch=compute_90,code=sm_90 \
    -c material_reference_cache.cu -o "$dest/reference_cache_gpu.o" > "$dest/build_cuda.log" 2>&1
mpicxx -O2 -std=c++17 -fopenmp -fno-fast-math -ffp-contract=off \
    -DXYZDBL -DGOTPM -DUSE_CUDA -I/opt/ohpc/pub/cuda/13.0.2/include \
    -c run_cell_reference_continuation_20261022.cpp -o "$dest/continuation.o" > "$dest/build_host.log" 2>&1
mpicxx -fopenmp "$dest/continuation.o" "$dest/reference_cache_gpu.o" \
    -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -o "$dest/continuation" > "$dest/build_link.log" 2>&1
sha256sum "$dest/continuation" "$dest/continuation.o" "$dest/reference_cache_gpu.o" > "$dest/binary.sha256"
sha256sum -c "$dest/source.sha256" > "$dest/postbuild_source_integrity.txt"
printf 'BUILD_PASS_NATIVE_NOT_EXECUTED_COMMON_STEP6_CLAIM_REQUIRED\n' > "$dest/status.txt"
