#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
dest=$root/material_compact_actual_source_v1
original=/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src
initial=$root/material_full_initial_v2_h200_mpi4/run
producer=$root/material_full_positive_rule_v1_ada_mpi4
transition=$root/material_full_positive_transition_v2_ada_mpi4
whole=$root/material_positive_rule_whole_audit_v1_mpi4
reference=$root/material_full_reference_v2_h200_mpi4/run/material_full_reference.json
accepted=$transition/run/continuous_step5/accepted_r
test "$#" = 3
first=$1
last=$2
chunk=$3
[[ "$first" =~ ^[0-9]+$ && "$last" =~ ^[0-9]+$ ]]
(( last > first && last - first <= 50000 ))
test "$(hostname)" = LagEunha
test "$(cat "$chunk/native_exit_code.txt")" = 0
canonical=$dest/power_edge_reconcile_${first}_${last}_v2
test "$(cat "$canonical/native_exit_code.txt")" = 0
run=$dest/power_edge_triangle_route_chunk_${first}_${last}_v1
source "$transition/run_env.sh"
export I_MPI_CXX=g++ I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0
export OMP_NUM_THREADS=4 OMP_PROC_BIND=false CUDA_VISIBLE_DEVICES=0
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY
cd "$dest"
sha256sum -c "$chunk/source_edge_wires.sha256" >/dev/null
sha256sum -c "$canonical/canonical.sha256" >/dev/null
sha256sum -c "$dest/weight_direction_v2/weights.sha256" >/dev/null
mkdir "$run"
sha256sum probe_power_triangle_route_v1.cpp run_power_triangle_route_lageunha_v1.sh \
  "$root/material_compact_rezone_parallel_v2.cpp" \
  "$root/material_compact_overlap_boundary_v2.h" \
  "$root/material_compact_overlap_cuts.h" \
  "$root/material_compact_overlap_triangle_route.h" \
  "$root/material_compact_overlap_events.h" \
  "$transition/Exam/libexam.a" \
  "$dest/weight_direction_v2/weights.csv" \
  "$chunk"/source_edge_manifest_rank{0,1,2,3}.csv \
  "$canonical"/canonical_owner{0,1,2,3}.csv > "$run/source.sha256"
mpicxx -std=gnu++17 -fext-numeric-literals -O2 -ffunction-sections -fdata-sections \
  -fopenmp -fno-fast-math -ffp-contract=off -DXYZDBL -DGOTPM -DUSE_CUDA \
  -I "$original/Exam" -I /opt/ohpc/pub/cuda/13.0.2/include \
  probe_power_triangle_route_v1.cpp -Wl,--gc-sections -L "$transition/Exam" -lexam \
  -L /opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lquadmath \
  -o "$run/probe_power_triangle_route"
sha256sum "$run/probe_power_triangle_route" > "$run/binary.sha256"
TIMEFORMAT='elapsed_seconds=%3R'
set +e
time timeout --signal=TERM --kill-after=60s 43200s \
  mpirun -np 4 "$run/probe_power_triangle_route" "$initial" \
  "$transition/immutable_payload.sha256" "$producer" "$whole" \
  "$reference" "$accepted" "$dest/context_v1" \
  "$dest/weight_direction_v2/weights.csv" 1 boundary_triangle_route_range \
  "$first" "$last" "$chunk" \
  </dev/null > "$run/native.log" 2>&1
result=$?
set -e
printf '%s\n' "$result" > "$run/native_exit_code.txt"
sha256sum -c "$run/source.sha256" > "$run/final_source_check.txt"
tail -25 "$run/native.log"
exit "$result"
