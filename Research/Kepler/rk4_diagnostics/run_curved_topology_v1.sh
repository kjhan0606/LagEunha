#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
dest=$root/rk4_diagnostics/curved_topology_$version
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
mkdir "$dest"
source "$root/run_env.sh"
export OMP_NUM_THREADS=2 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
for h in gfs_material_map12.h gfs_material_p3.h gfs_material_p3_radial_input.h gfs_material_triangle.h gfs_moment.h; do cp "$root/rk4_diagnostics/material_p3_finite_v2/$h" "$dest/"; done
for f in material_curved_packet.h material_curved_bounds.h material_curved_compiled.h material_curved_roots.h material_curved_trace_bounds.h material_curved_partition.h material_curved_topology.h test_curved_transfer_v1.cpp test_curved_topology_v1.cpp; do cp "$root/rk4_diagnostics/$f" "$dest/"; done
packet=$root/rk4_diagnostics/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
sites=$root/rk4_diagnostics/curved_map_algebra_v2/advected_sites.csv
sha256sum "$packet" "$sites" > "$dest/inputs.sha256"
mpicxx -O2 -std=c++17 -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_curved_topology_v1.cpp" -o "$dest/contract.exe" > "$dest/build.log" 2>&1
mpirun -np 4 "$dest/contract.exe" "$packet" "$sites" "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
printf 'AUDIT_COMPLETE\n' > "$dest/status.txt"
