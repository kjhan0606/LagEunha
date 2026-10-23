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
run=$dest/context_v1
test "$(hostname)" = LagEunha
source "$transition/run_env.sh"
export I_MPI_CXX=g++ I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0
export OMP_NUM_THREADS=4 OMP_PROC_BIND=false CUDA_VISIBLE_DEVICES=0
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY
cd "$dest"
test -s lageunha_immutable_input_check.txt
test -s lageunha_positive_input_check.txt
sha256sum -c lageunha_actual_source_input.sha256 > context_previous_input_check.txt
mkdir "$run"
sha256sum probe_context.cpp "$transition/Exam/libexam.a" \
  "$reference" "$whole/contract.json" \
  "$root/material_full_reference_v2_h200_mpi4/run/material_native_geometry_r0.bin" \
  "$root/material_full_reference_v2_h200_mpi4/run/material_native_geometry_r1.bin" \
  "$root/material_full_reference_v2_h200_mpi4/run/material_native_geometry_r2.bin" \
  "$root/material_full_reference_v2_h200_mpi4/run/material_native_geometry_r3.bin" \
  > context_source_input.sha256
mpicxx -std=gnu++17 -fext-numeric-literals -O2 -ffunction-sections -fdata-sections \
  -fopenmp -fno-fast-math -ffp-contract=off -DXYZDBL -DGOTPM -DUSE_CUDA \
  -I "$original/Exam" -I /opt/ohpc/pub/cuda/13.0.2/include \
  probe_context.cpp -Wl,--gc-sections -L "$transition/Exam" -lexam \
  -L /opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lquadmath -o probe_context_lageunha
sha256sum probe_context_lageunha > context_binary.sha256
TIMEFORMAT='elapsed_seconds=%3R'
set +e
time timeout --signal=TERM --kill-after=60s 7200s \
  mpirun -np 4 ./probe_context_lageunha "$initial" \
  "$transition/immutable_payload.sha256" "$producer" "$whole" \
  "$reference" "$accepted" "$run" </dev/null > context_native.log 2>&1
result=$?
set -e
printf '%s\n' "$result" > context_native_exit_code.txt
sha256sum -c context_source_input.sha256 > context_final_source_input_check.txt
if test -f "$run/geometry_r0.bin"; then
  sha256sum "$run"/geometry_r*.bin > context_geometry.sha256
fi
tail -35 context_native.log
exit "$result"
