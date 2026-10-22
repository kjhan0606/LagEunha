#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
guard=${2:-$root/rk4_diagnostics/curved_events_v9}
ranks=${3:-4}
triangle=${4:-}
dest=$root/rk4_diagnostics/curved_complement_$version
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
test "$(cat "$guard/status.txt")" = AUDIT_COMPLETE
mkdir "$dest"
source "${GFS_CURVED_FROZEN_ENV:-$root/run_env.sh}"
export OMP_NUM_THREADS=2 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
sha256sum -c "$guard/inputs.sha256" > "$dest/input_integrity.txt"
for f in gfs_material_map12.h gfs_material_p3.h gfs_material_p3_radial_input.h gfs_material_triangle.h gfs_moment.h material_curved_packet.h material_curved_bounds.h material_curved_compiled.h material_curved_roots.h material_curved_trace_bounds.h material_curved_partition.h material_curved_topology.h material_curved_events.h test_curved_transfer_v1.cpp; do cp "$guard/$f" "$dest/"; done
cp "$root/rk4_diagnostics/test_curved_complement_v1.cpp" "$dest/"
# Reuse the same compiled map/chart polynomials and prune only targets excluded
# by a whole-source interval certificate. Field/root/bound headers stay frozen.
cp "$root/rk4_diagnostics/material_curved_partition.h" "$dest/"
cp "$root/rk4_diagnostics/test_curved_transfer_v1.cpp" "$dest/"
sha256sum "$dest/"*.h "$dest/"*.cpp > "$dest/source.sha256"
mpicxx -O2 -std=c++17 -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_curved_complement_v1.cpp" -o "$dest/contract.exe" > "$dest/build.log" 2>&1
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
printf '%s\n' "$guard" > "$dest/event_guard.txt"
packet=$root/rk4_diagnostics/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
sites=$root/rk4_diagnostics/curved_map_algebra_v2/advected_sites.csv
set +e
args=("$packet" "$sites" "$guard" "$dest")
if test -n "$triangle"; then args+=("$triangle"); fi
mpirun -np "$ranks" "$dest/contract.exe" "${args[@]}" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
if test "$result" = 0; then printf 'ESTIMATE_AUDIT_COMPLETE_NOT_COLD_ACCEPTANCE\n' > "$dest/status.txt"; else printf 'FAILED_EXIT_%s\n' "$result" > "$dest/status.txt"; fi
exit "$result"
