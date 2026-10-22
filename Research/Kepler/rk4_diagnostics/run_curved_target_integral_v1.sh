#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
case_id=${2:-225}
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
case "$case_id" in 225|137) ;; *) exit 2;; esac
dest=$root/rk4_diagnostics/curved_target_integral_${version}_${case_id}
guard=$root/rk4_diagnostics/curved_events_v9
mkdir "$dest"
source "${GFS_CURVED_FROZEN_ENV:-$root/run_env.sh}"
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
for f in gfs_material_map12.h gfs_material_p3.h gfs_material_p3_radial_input.h gfs_material_triangle.h gfs_moment.h; do cp "$guard/$f" "$dest/"; done
for f in material_curved_packet.h material_curved_bounds.h material_curved_compiled.h material_curved_roots.h material_curved_trace_bounds.h material_curved_partition.h material_curved_topology.h material_curved_events.h material_curved_branches.h material_curved_branch_bounds.h material_forward_ball.h material_original_oracle.h material_native_forward.h material_target_inverse.h material_target_general.h material_target_domain.h test_curved_transfer_v1.cpp test_curved_complement_v1.cpp test_curved_branch_v1.cpp test_curved_original_oracle_v1.cpp test_curved_target_integral_v1.cpp; do cp "$root/rk4_diagnostics/$f" "$dest/"; done
sha256sum "$dest/"*.h "$dest/"*.cpp > "$dest/source.sha256"
sha256sum -c "$guard/inputs.sha256" > "$dest/input_integrity.txt"
sha256sum "$root/rk4_diagnostics/curved_fold_transfer_v5_${case_id}/positive_points.bin" "$root/rk4_diagnostics/curved_fold_transfer_v5_${case_id}/target_contribution_bounds.csv" > "$dest/retained_inputs.sha256"
printf '%s\n' "${GFS_CURVED_LAUNCH_ENVELOPE:-missing}" > "$dest/runner_envelope.txt"
set +e
mpicxx -O2 -std=c++17 -fext-numeric-literals -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_curved_target_integral_v1.cpp" -lquadmath -o "$dest/contract.exe" > "$dest/build.log" 2>&1
build_result=$?
set -e
printf '%s\n' "$build_result" > "$dest/build_exit_code.txt"
if test "$build_result" != 0; then printf 'BUILD_FAILED_EXIT_%s\n' "$build_result" > "$dest/status.txt"; exit "$build_result"; fi
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
set +e
mpirun -np 1 "$dest/contract.exe" "$root/rk4_diagnostics/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin" "$root/rk4_diagnostics/curved_map_algebra_v2/advected_sites.csv" "$guard" "$dest" "$case_id" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/retained_inputs.sha256" > "$dest/final_retained_input_integrity.txt"
if test -f "$dest/physical_positive_points.bin"; then sha256sum "$dest/physical_positive_points.bin" > "$dest/new_positive_packet.sha256"; fi
if test "$result" = 0; then printf 'PHYSICAL_TARGET_QUADRATURE_COMPLETE_ERROR_LEDGER_PENDING\n' > "$dest/status.txt"; else printf 'FAILED_EXIT_%s\n' "$result" > "$dest/status.txt"; fi
exit "$result"
