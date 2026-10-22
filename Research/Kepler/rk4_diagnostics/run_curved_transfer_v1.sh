#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
dest=$root/rk4_diagnostics/curved_transfer_$version
case "$(hostname)" in LagEunha|lageunha) ;; *) printf 'Use authorized lageunha\n' >&2;exit 2;; esac
mkdir "$dest"
source "$root/run_env.sh"
export OMP_NUM_THREADS=2 I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
export OMP_PROC_BIND=false
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY
unset I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
for h in gfs_material_map12.h gfs_material_p3.h gfs_material_p3_radial_input.h gfs_material_triangle.h gfs_moment.h; do cp "$root/rk4_diagnostics/material_p3_finite_v2/$h" "$dest/"; done
for f in material_curved_packet.h material_curved_bounds.h material_curved_compiled.h material_curved_roots.h material_curved_trace_bounds.h material_curved_partition.h test_curved_transfer_v1.cpp; do cp "$root/rk4_diagnostics/$f" "$dest/"; done
packet=$root/rk4_diagnostics/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
sites=$root/rk4_diagnostics/curved_map_algebra_v2/advected_sites.csv
sha256sum "$packet" "$sites" > "$dest/inputs.sha256"
mpicxx -O2 -std=c++17 -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_curved_transfer_v1.cpp" -o "$dest/contract.exe" > "$dest/build.log" 2>&1
for tri in 182 149; do
    mkdir "$dest/debug_triangle_$tri"
    printf '%s\n' "$tri" > "$dest/debug_triangle_$tri/debug_filter.txt"
    mpirun -np 1 "$dest/contract.exe" "$packet" "$sites" "$dest/debug_triangle_$tri" 12 10 0 0 "$tri" > "$dest/debug_triangle_$tri/contract.json" 2> "$dest/debug_triangle_$tri/stderr.log"
done
for mode in coarse_mpi1 coarse_mpi4 fine_mpi4 weighted_mpi4; do
    ranks=4;nr=12;na=10;weighted=0;save=0
    if test "$mode" = coarse_mpi1; then ranks=1;fi
    if test "$mode" = fine_mpi4; then nr=24;na=16;save=1;fi
    if test "$mode" = weighted_mpi4; then nr=24;na=16;weighted=1;fi
    mkdir "$dest/$mode"
    mpirun -np "$ranks" "$dest/contract.exe" "$packet" "$sites" "$dest/$mode" "$nr" "$na" "$weighted" "$save" > "$dest/$mode/contract.json" 2> "$dest/$mode/stderr.log"
done
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
cmp "$dest/coarse_mpi1/transfer_moments.csv" "$dest/coarse_mpi4/transfer_moments.csv"
cmp "$dest/coarse_mpi1/angular_uncertainty.csv" "$dest/coarse_mpi4/angular_uncertainty.csv"
printf 'AUDIT_COMPLETE\n' > "$dest/status.txt"
