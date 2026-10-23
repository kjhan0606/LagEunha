#!/bin/bash
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
run=$dest/power_v3
test "$(hostname)" = LagEunha
source "$transition/run_env.sh"
export I_MPI_CXX=g++ I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0
export OMP_NUM_THREADS=8 OMP_PROC_BIND=false CUDA_VISIBLE_DEVICES=0
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY
cd "$dest"
test "$(cat context_native_exit_code.txt)" = 0
sha256sum -c context_geometry.sha256 > power_cached_geometry_check.txt
sha256sum -c context_source_input.sha256 > power_previous_context_input_check.txt
sha256sum -c lageunha_actual_source_input.sha256 > power_accepted_input_check.txt
mkdir "$run"
sha256sum probe_power.cpp "$root/material_compact_rezone.cpp" \
  "$root/material_compact_chart_topology.h" \
  "$transition/Exam/libexam.a" > "$run/source.sha256"
mpicxx -std=gnu++17 -fext-numeric-literals -O2 -ffunction-sections -fdata-sections \
  -fopenmp -fno-fast-math -ffp-contract=off -DXYZDBL -DGOTPM -DUSE_CUDA \
  -I "$original/Exam" -I /opt/ohpc/pub/cuda/13.0.2/include \
  probe_power.cpp -Wl,--gc-sections -L "$transition/Exam" -lexam \
  -L /opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lquadmath \
  -o "$run/probe_power"
sha256sum "$run/probe_power" > "$run/binary.sha256"
TIMEFORMAT='elapsed_seconds=%3R'
set +e
time timeout --signal=TERM --kill-after=60s 7200s \
  mpirun -np 4 "$run/probe_power" "$initial" \
  "$transition/immutable_payload.sha256" "$producer" "$whole" \
  "$reference" "$accepted" "$dest/context_v1" \
  </dev/null > "$run/native.log" 2>&1
result=$?
set -e
printf '%s\n' "$result" > "$run/native_exit_code.txt"
sha256sum -c "$run/source.sha256" > "$run/final_source_check.txt"
tail -35 "$run/native.log"
exit "$result"
