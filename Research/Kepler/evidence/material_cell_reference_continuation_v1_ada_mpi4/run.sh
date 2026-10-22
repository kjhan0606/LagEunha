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
printf 'CHECKING_SOURCE_AND_INPUTS_BEFORE_REFERENCE_CONTINUATION\n' > "$dest/status.txt"
sha256sum -c "$dest/source.sha256" > "$dest/prelaunch_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/prelaunch_input_integrity.txt"
sha256sum -c "$dest/binary.sha256" > "$dest/prelaunch_binary_integrity.txt"
active=$(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader)
if test -n "$active"; then printf '%s\n' "$active" > "$dest/device_busy.txt"; exit 3; fi
python3 - "$dest" <<'CHECK'
from pathlib import Path
import sys,shutil
memory=dict((x.split(':',1)[0],x.split(':',1)[1].strip()) for x in Path('/proc/meminfo').read_text().splitlines())
assert int(memory['MemAvailable'].split()[0])*1024>128*2**30
assert shutil.disk_usage(sys.argv[1]).free>8*2**30
(Path(sys.argv[1])/'host_memory_at_launch.txt').write_text(Path('/proc/meminfo').read_text())
CHECK
nvidia-smi --query-gpu=index,uuid,name,memory.total,memory.used --format=csv > "$dest/device.txt"
mkdir "$dest/run"
printf 'RUNNING_REFERENCE_CACHE_AND_TWO_CONTINUOUS_NATIVE_STEPS\n' > "$dest/status.txt"
set +e
timeout --signal=TERM --kill-after=30s 1800s mpirun -np 4 "$dest/continuation" "$root/material_full_initial_v2_h200_mpi4/run" "$root/material_full_positive_rule_v1_ada_mpi4" "$root/material_full_positive_transition_v2_ada_mpi4/immutable_payload.sha256" "$root/material_full_positive_transition_v2_ada_mpi4/run/continuous_step5/accepted_r" "$root/material_cell_condensation_structure_v1_mpi4/run" "$root/material_cell_condensation_factor_v1_mpi4/run" "$root/material_cell_condensation_actual_rhs_v1_ada_mpi4/run" "$dest/reference_source_binding.txt" "$dest/run" > "$dest/native.log" 2>&1
native_result=$?
set -e
printf '%s\n' "$native_result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/final_input_integrity.txt"
sha256sum -c "$dest/binary.sha256" > "$dest/final_binary_integrity.txt"
if test "$native_result" = 0; then
 printf 'TWO_NATIVE_STEPS_COMPLETE_LONG_RELEASE_NOT_AUTHORIZED_BY_GATES\n' > "$dest/status.txt"
else
 printf 'REFERENCE_CACHE_CONTINUATION_REJECTED_%s\n' "$native_result" > "$dest/status.txt"
fi
exit "$native_result"
