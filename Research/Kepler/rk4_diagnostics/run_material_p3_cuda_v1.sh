#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
dest=$root/material_p3_cuda_v1
case "$(hostname)" in LagEunha|lageunha) ;; *) printf 'Use authorized lageunha\n' >&2;exit 2;; esac
active=$(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader)
if test -n "$active"; then printf 'GPU compute process present; no launch\n' >&2;exit 3;fi
mkdir "$dest"
source "$root/../run_env.sh"
export CUDA_VISIBLE_DEVICES=0
nvidia-smi --query-gpu=index,name,memory.total,memory.used --format=csv > "$dest/device.txt"
cp "$root/material_p3_hybrid_v2/"gfs_material_p3{,_radial_input}.h "$root/material_p3_hybrid_v2/gfs_material_triangle.h" "$dest/"
cp "$root/material_p3_hybrid_v2/mpi1/level0_remainder.csv" "$root/material_p3_hybrid_v2/mpi1/level1_remainder.csv" "$root/test_material_p3_cuda_v1.cu" "$dest/"
nvcc -O2 -std=c++17 -I"$dest" "$dest/test_material_p3_cuda_v1.cu" -o "$dest/contract.exe" > "$dest/build.log" 2>&1
"$dest/contract.exe" "$dest/level0_remainder.csv" "$dest/level1_remainder.csv" > "$dest/contract.json" 2> "$dest/stderr.log"
sha256sum "$dest/contract.exe" "$dest/test_material_p3_cuda_v1.cu" "$dest/gfs_material_p3.h" "$dest/gfs_material_p3_radial_input.h" "$dest/gfs_material_triangle.h" > "$dest/artifact.sha256"
printf 'AUDIT_COMPLETE\n' > "$dest/status.txt"
