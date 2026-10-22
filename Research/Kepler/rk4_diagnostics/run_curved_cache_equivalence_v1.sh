#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
dest=$root/rk4_diagnostics/curved_cache_equivalence_v1
guard=$root/rk4_diagnostics/curved_events_v9
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
mkdir "$dest"
source "${GFS_CURVED_FROZEN_ENV:-$root/run_env.sh}"
export OMP_NUM_THREADS=2 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
for mode in old new; do
    if test "$mode" = old; then parent=$root/rk4_diagnostics/curved_complement_origin_v1; else parent=$root/rk4_diagnostics/curved_complement_cache_origin_v2; fi
    out=$dest/$mode
    mkdir "$out"
    for f in gfs_material_map12.h gfs_material_p3.h gfs_material_p3_radial_input.h gfs_material_triangle.h gfs_moment.h material_curved_packet.h material_curved_bounds.h material_curved_compiled.h material_curved_roots.h material_curved_trace_bounds.h material_curved_partition.h material_curved_topology.h material_curved_events.h test_curved_transfer_v1.cpp test_curved_complement_v1.cpp; do cp "$parent/$f" "$out/"; done
    cp "$root/rk4_diagnostics/test_curved_cache_equivalence_v1.cpp" "$out/"
    sha256sum "$out/"*.h "$out/"*.cpp > "$out/source.sha256"
    mpicxx -O2 -std=c++17 -fopenmp -fno-fast-math -ffp-contract=off -I"$out" "$out/test_curved_cache_equivalence_v1.cpp" -o "$out/contract.exe" > "$out/build.log" 2>&1
    mpirun -np 1 "$out/contract.exe" "$root/rk4_diagnostics/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin" "$root/rk4_diagnostics/curved_map_algebra_v2/advected_sites.csv" "$guard" "$out" 133 > "$out/contract.json" 2> "$out/stderr.log"
done
cmp "$dest/old/moments.csv" "$dest/new/moments.csv" > "$dest/moments_cmp.txt"
cmp "$dest/old/points.bin" "$dest/new/points.bin" > "$dest/points_cmp.txt"
printf 'AUDIT_COMPLETE\n' > "$dest/status.txt"
