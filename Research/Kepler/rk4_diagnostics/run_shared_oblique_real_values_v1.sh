#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
diag=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
dest=$diag/shared_oblique_real_values_${1:-v1}_mpi1
mkdir "$dest"
source "$GFS_CURVED_FROZEN_ENV"
export OMP_NUM_THREADS=16 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
previous=$diag/shared_real_graph_values_v1_mpi1
for input in "$previous/"*.h "$previous/"*.cpp; do cp "$input" "$dest/"; done
cp "$diag/test_shared_oblique_real_values_v1.cpp" "$dest/"
sha256sum -c "$previous/source.sha256" > "$dest/previous_source_integrity.txt"
sha256sum -c "$previous/output.sha256" > "$dest/previous_output_integrity.txt"
inputs=$diag/shared_mixed_runtime_inputs_v1
sha256sum -c "$inputs/input.sha256" > "$dest/registry_input_integrity.txt"
sha256sum -c "$inputs/output.sha256" > "$dest/registry_output_integrity.txt"
packet=$diag/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
saved=$diag/shared_trace_rhs_v1_mpi4
sites=$diag/curved_map_algebra_v2/advected_sites.csv
weights=$diag/curved_events_v9/site_weights.csv
slabs=$diag/shared_oblique_slabs_v1_mpi1/actual_oblique_slabs_rank0.bin
pressure=$diag/geometric_pressure_v2_mpi4
sha256sum "$dest/"*.h "$dest/"*.cpp > "$dest/source.sha256"
sha256sum "$packet" "$saved/shared_state.bin" "$sites" "$weights" "$slabs" "$diag/shared_atlas_adoption_v7/typed_packet_registry.csv" "$inputs/real_graph_groups.csv" "$inputs/panel_jobs.csv" "$pressure/geometric_cuda.o" "$diag/SHARED_REAL_GRAPH_VALUES_POLICY_20261022.md" > "$dest/input.sha256"
cat "$inputs/input.sha256" "$inputs/output.sha256" >> "$dest/input.sha256"
printf '%s\n' "$GFS_CURVED_LAUNCH_ENVELOPE" > "$dest/runner_envelope.txt"
set +e
mpicxx -O2 -std=c++17 -fext-numeric-literals -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_shared_oblique_real_values_v1.cpp" "$pressure/geometric_cuda.o" -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lquadmath -o "$dest/contract.exe" > "$dest/build.log" 2>&1
build=$?
set -e
printf '%s\n' "$build" > "$dest/build_exit_code.txt"
if test "$build" != 0; then printf 'BUILD_FAILED\n' > "$dest/status.txt"; exit 1; fi
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
set +e
mpirun -np 1 "$dest/contract.exe" "$packet" "$saved" "$sites" "$weights" "$inputs" "$slabs" "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/final_input_integrity.txt"
shopt -s nullglob
outputs=("$dest/"*.json "$dest/"*.csv "$dest/"*.bin)
if test "${#outputs[@]}" -gt 0; then sha256sum "${outputs[@]}" > "$dest/output.sha256"; fi
printf 'OBLIQUE_REAL_VALUES_EXIT_%s\n' "$result" > "$dest/status.txt"
exit "$result"
