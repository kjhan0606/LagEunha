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
test "$(hostname)" = LagEunha
source "$transition/run_env.sh"
export I_MPI_CXX=g++ I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0
export OMP_NUM_THREADS=2 OMP_PROC_BIND=false
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY
cd "$dest"
sha256sum -c "$transition/immutable_payload.sha256" > lageunha_immutable_input_check.txt
sha256sum -c "$transition/positive_rule_input.sha256" > lageunha_positive_input_check.txt
sha256sum "$dest/probe.cpp" "$root/material_compact_rezone.cpp" \
  "$reference" "$whole/contract.json" \
  "$transition/run/continuous_step5/step.json" \
  "${accepted}0.bin" "${accepted}1.bin" "${accepted}2.bin" "${accepted}3.bin" \
  > lageunha_actual_source_input.sha256
mpicxx -std=gnu++17 -fext-numeric-literals -O2 -ffunction-sections -fdata-sections \
  -fopenmp -fno-fast-math -ffp-contract=off -DXYZDBL -DGOTPM -DUSE_CUDA \
  -I "$original/Exam" -I /opt/ohpc/pub/cuda/13.0.2/include \
  probe.cpp -Wl,--gc-sections -L /opt/ohpc/pub/cuda/13.0.2/lib64 \
  -lquadmath -lcudart -o probe_lageunha
sha256sum probe_lageunha > lageunha_probe_binary.sha256
TIMEFORMAT='elapsed_seconds=%3R'
set +e
time timeout --signal=TERM --kill-after=60s 7200s \
  mpirun -np 4 ./probe_lageunha "$initial" \
  "$transition/immutable_payload.sha256" "$producer" "$whole" \
  "$reference" "$accepted" </dev/null > lageunha_native.log 2>&1
result=$?
set -e
printf '%s\n' "$result" > lageunha_native_exit_code.txt
sha256sum -c lageunha_actual_source_input.sha256 > lageunha_final_source_input_check.txt
tail -35 lageunha_native.log
exit "$result"
