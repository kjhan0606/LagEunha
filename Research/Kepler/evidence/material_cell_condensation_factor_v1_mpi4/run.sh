#!/usr/bin/env bash
set -euo pipefail
dest=$(cd -- "$(dirname -- "$0")" && pwd)
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
test "$(hostname)" = LagEunha
test ! -e "$dest/run"
source "$dest/run_env.sh"
export I_MPI_CXX=g++ I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 OMP_NUM_THREADS=8 OMP_PROC_BIND=false
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
sha256sum -c "$dest/source_revision2.sha256" > "$dest/prelaunch_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/prelaunch_input_integrity.txt"
sha256sum -c "$dest/binary_revision2.sha256" > "$dest/prelaunch_binary_integrity.txt"
mkdir "$dest/run"
printf 'RUNNING_CPU_ORIGINAL_CELL_FACTORS_NO_EVOLUTION\n' > "$dest/status.txt"
set +e
timeout --signal=TERM --kill-after=30s 28800s mpirun -np 4 "$dest/build_cell_condensation_revision2" "$root/material_full_initial_v2_h200_mpi4/run" "$root/material_cell_condensation_structure_v1_mpi4/run" "$dest/setup_capacity.json" "$dest/run" > "$dest/native.log" 2>&1
native_result=$?
set -e
printf '%s\n' "$native_result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source_revision2.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/final_input_integrity.txt"
sha256sum -c "$dest/binary_revision2.sha256" > "$dest/final_binary_integrity.txt"
if test "$native_result" = 0; then
 printf 'ORIGINAL_CELL_FACTORS_COMPLETE_ACTUAL_RHS_GPU_PENDING\n' > "$dest/status.txt"
else
 printf 'ORIGINAL_CELL_FACTORS_REJECTED_%s\n' "$native_result" > "$dest/status.txt"
fi
exit "$native_result"
