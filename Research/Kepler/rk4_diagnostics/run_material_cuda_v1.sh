#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
src=/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src
dest=$root/rk4_diagnostics/material_cuda_v1
case "$(hostname)" in LagEunha|lageunha) ;; *) printf 'Use authorized lageunha\n' >&2;exit 2;; esac
if [ -n "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ]; then printf 'GPU occupied; audit not launched\n' >&2;exit 3;fi
mkdir "$dest"
source "$root/run_env.sh"
export OMP_NUM_THREADS=1
cd "$src"
sha256sum Exam/gfs_material_map.h > "$dest/source.sha256"
cp Exam/gfs_material_map.h "$dest/"
cp "$root/rk4_diagnostics/test_material_cuda_v1.cu" "$dest/"
nvcc -O2 -std=c++17 -arch=sm_89 -I"$dest" "$dest/test_material_cuda_v1.cu" -o "$dest/contract.exe" > "$dest/build.log" 2>&1
"$dest/contract.exe" > "$dest/contract.json" 2> "$dest/stderr.log"
sha256sum -c "$dest/source.sha256" > "$dest/source_unchanged.txt"
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
printf 'AUDIT_COMPLETE\n' > "$dest/status.txt"
