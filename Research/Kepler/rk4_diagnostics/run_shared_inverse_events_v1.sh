#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
diag=$root/rk4_diagnostics
version=${1:-v1}
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
dest=$diag/shared_inverse_events_${version}_mpi1
mkdir "$dest"
source "$GFS_CURVED_FROZEN_ENV"
export OMP_NUM_THREADS=4 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
witness=$diag/shared_atlas_witness_v1_mpi4
repair=$diag/shared_atlas_repairs_v1_mpi1
budget_policy=$diag/SHARED_ATLAS_BUDGET_POLICY_20261021.json
budget=$(python3 -c 'import json,sys; p=json.load(open(sys.argv[1])); assert p["increase_factor"]==16 and p["fixed_gauss_order_per_axis"]==16 and not p["full_quadrature_allowed"]; print(p["new_per_target_query_limit"])' "$budget_policy")
python3 -c 'import json,sys; p=json.load(open(sys.argv[1])); assert p["limited_saved_repair_audit_pass"] and p["shared_node_inverse_interior_and_contraction_inequalities_pass"]' "$diag/SHARED_ATLAS_REPAIRS_20261021.json"
for input in "$repair/"*.h "$repair/"*.cpp; do cp "$input" "$dest/"; done
cp "$diag/test_shared_inverse_events_v1.cpp" "$dest/"
sha256sum -c "$witness/source.sha256" > "$dest/original_witness_source_integrity.txt"
sha256sum -c "$witness/output.sha256" > "$dest/original_witness_output_integrity.txt"
sha256sum -c "$repair/output.sha256" > "$dest/certified_inverse_integrity.txt"
packet=$diag/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
saved=$diag/shared_trace_rhs_v1_mpi4
sites=$diag/curved_map_algebra_v2/advected_sites.csv
weights=$diag/curved_events_v9/site_weights.csv
actual=$diag/shared_registry_domains_v1_mpi4
pressure=$diag/geometric_pressure_v2_mpi4
sha256sum "$dest/"*.h "$dest/"*.cpp > "$dest/source.sha256"
sha256sum "$packet" "$saved/shared_state.bin" "$sites" "$weights" "$witness/failed_charts_rank"*.bin "$witness/chart_events_rank"*.bin "$repair/vertex_seed_inverse_rank0.bin" "$actual/domain_nodes_rank"*.csv "$pressure/geometric_cuda.o" "$budget_policy" "$diag/ASTRA_SHARED_FAILED_ATLAS_END_REVIEW_20261021.md" "$diag/ASTRA_SHARED_FAILED_ATLAS_BUDGET_ADDENDUM_20261021.md" > "$dest/input.sha256"
printf '%s\n' "$GFS_CURVED_LAUNCH_ENVELOPE" > "$dest/runner_envelope.txt"
set +e
mpicxx -O2 -std=c++17 -fext-numeric-literals -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_shared_inverse_events_v1.cpp" "$pressure/geometric_cuda.o" -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lquadmath -o "$dest/contract.exe" > "$dest/build.log" 2>&1
build=$?
set -e
printf '%s\n' "$build" > "$dest/build_exit_code.txt"
if test "$build" != 0; then printf 'BUILD_FAILED\n' > "$dest/status.txt"; exit 1; fi
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
set +e
mpirun -np 1 "$dest/contract.exe" "$packet" "$saved" "$sites" "$weights" "$witness" "$repair" "$actual" "$dest" "$budget" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/final_input_integrity.txt"
sha256sum "$dest/"*.json "$dest/"*.csv "$dest/"*.bin > "$dest/output.sha256"
printf 'SHARED_INVERSE_EVENT_BINDING_EXIT_%s\n' "$result" > "$dest/status.txt"
exit "$result"
