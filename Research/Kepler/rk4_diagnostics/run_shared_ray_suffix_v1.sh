#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
diag=$root/rk4_diagnostics
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
dest=$diag/shared_ray_suffix_${1:-v1}_mpi1
mkdir "$dest"
source "$GFS_CURVED_FROZEN_ENV"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
old=$diag/shared_origin_atlas_v2_mpi1
arrangements=$diag/shared_inverse_events_v1_mpi1
python3 -c 'import json,sys; p=json.load(open(sys.argv[1])); assert p["both_original_failure_disks_in_positive_ray_square_half_plane_after_collection"]' "$diag/SHARED_RAY_DEPENDENCY_20261021.json"
for input in "$old/"*.h "$old/"*.cpp; do cp "$input" "$dest/"; done
cp "$diag/test_shared_ray_suffix_v1.cpp" "$diag/material_shared_collected_radius.h" "$dest/"
sha256sum -c "$old/source.sha256" > "$dest/retained_source_integrity.txt"
sha256sum -c "$old/output.sha256" > "$dest/retained_output_integrity.txt"
packet=$diag/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
saved=$diag/shared_trace_rhs_v1_mpi4
sites=$diag/curved_map_algebra_v2/advected_sites.csv
weights=$diag/curved_events_v9/site_weights.csv
policy=$diag/SHARED_ATLAS_BUDGET_POLICY_20261021.json
pressure=$diag/geometric_pressure_v2_mpi4
cap=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["new_per_target_query_limit"])' "$policy")
sha256sum "$dest/"*.h "$dest/"*.cpp > "$dest/source.sha256"
sha256sum "$packet" "$saved/shared_state.bin" "$sites" "$weights" "$old/graded_source133_atlas_rank0.bin" "$old/origin_atlas_cost_and_clearance.csv" "$arrangements/new_failed_chart_arrangements_rank0.bin" "$pressure/geometric_cuda.o" "$policy" "$diag/SHARED_RAY_DEPENDENCY_20261021.json" > "$dest/input.sha256"
printf '%s\n' "$GFS_CURVED_LAUNCH_ENVELOPE" > "$dest/runner_envelope.txt"
set +e
mpicxx -O2 -std=c++17 -fext-numeric-literals -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_shared_ray_suffix_v1.cpp" "$pressure/geometric_cuda.o" -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lquadmath -o "$dest/contract.exe" > "$dest/build.log" 2>&1
build=$?
set -e
printf '%s\n' "$build" > "$dest/build_exit_code.txt"
if test "$build" != 0; then printf 'BUILD_FAILED\n' > "$dest/status.txt"; exit 1; fi
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
set +e
mpirun -np 1 "$dest/contract.exe" "$packet" "$saved" "$sites" "$weights" "$old" "$arrangements" "$dest" "$cap" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/final_input_integrity.txt"
sha256sum "$dest/"*.json "$dest/"*.csv "$dest/"*.bin > "$dest/output.sha256"
printf 'COLLECTED_RAY_SUFFIX_EXIT_%s\n' "$result" > "$dest/status.txt"
exit "$result"
