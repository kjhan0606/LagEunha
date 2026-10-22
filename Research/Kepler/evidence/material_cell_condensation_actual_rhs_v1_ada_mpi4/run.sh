#!/usr/bin/env bash
set -euo pipefail
dest=$(cd -- "$(dirname -- "$0")" && pwd)
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
test "$(hostname)" = LagEunha
test ! -e "$dest/run"
source "$dest/run_env.sh"
export I_MPI_CXX=g++ I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 OMP_NUM_THREADS=8 OMP_PROC_BIND=false CUDA_VISIBLE_DEVICES=0
export GFS_MATERIAL_DEVICE_MASS=1 GFS_MATERIAL_SAVE_SOURCE_CONTEXT=1
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
sha256sum -c "$dest/source.sha256" > "$dest/prelaunch_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/prelaunch_input_integrity.txt"
sha256sum -c "$dest/binary.sha256" > "$dest/prelaunch_binary_integrity.txt"
active=$(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader)
if test -n "$active"; then printf '%s\n' "$active" > "$dest/device_busy.txt"; exit 3; fi
nvidia-smi --query-gpu=index,uuid,name,memory.total,memory.used --format=csv > "$dest/device.txt"
mkdir "$dest/run"
printf 'RUNNING_FIRST_ACTUAL_STEP6_RHS_NO_ACCEPTED_REPLAY\n' > "$dest/status.txt"
set +e
timeout --signal=TERM --kill-after=30s 1800s mpirun -np 4 "$dest/actual_rhs" "$root/material_full_initial_v2_h200_mpi4/run" "$root/material_full_positive_rule_v1_ada_mpi4" "$root/material_full_positive_transition_v2_ada_mpi4/immutable_payload.sha256" "$root/material_full_positive_transition_v2_ada_mpi4/run/continuous_step5/accepted_r" "$root/material_cell_condensation_structure_v1_mpi4/run" "$root/material_cell_condensation_factor_v1_mpi4/run" "$dest/run" > "$dest/native.log" 2>&1
native_result=$?
set -e
printf '%s\n' "$native_result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/final_input_integrity.txt"
sha256sum -c "$dest/binary.sha256" > "$dest/final_binary_integrity.txt"
if test "$native_result" = 0; then
 printf 'ACTUAL_NEXT_RHS_AND_CONDENSATION_PASS_RK4_CACHE_PENDING\n' > "$dest/status.txt"
else
 printf 'ACTUAL_NEXT_RHS_REJECTED_%s\n' "$native_result" > "$dest/status.txt"
fi
exit "$native_result"
