#!/usr/bin/env bash
#SBATCH --job-name=kepler-cell-reference
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=192G
#SBATCH --time=00:40:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/material_cell_reference_continuation_v1_h200_mpi4/slurm-%j.out
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
dest=$root/material_cell_reference_continuation_v1_h200_mpi4
donor=$root/material_cell_reference_continuation_v1_ada_mpi4
claim=$root/material_step6_reference_continuation.claim
native_result=NOT_STARTED
claim_owned=0
finish() {
    runner_result=$?
    trap - EXIT
    printf '%s\n' "$runner_result" > "$dest/runner_exit_code.txt"
    if test "$claim_owned" = 1; then
        printf 'job=%s\nrunner_exit=%s\nnative_exit=%s\n' "$SLURM_JOB_ID" "$runner_result" "$native_result" > "$claim/terminal.txt"
    fi
    exit "$runner_result"
}
trap finish EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
test ! -e "$dest/run"
source "$dest/run_env.sh"
export I_MPI_CXX=g++ I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0
export OMP_NUM_THREADS=8 OMP_PROC_BIND=false
export GFS_MATERIAL_DEVICE_MASS=1 GFS_MATERIAL_SAVE_SOURCE_CONTEXT=1
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
printf '%s\n' "$SLURM_JOB_ID" > "$dest/slurm_job_id.txt"
printf 'CHECKING_SOURCE_AND_INPUTS_BEFORE_H200_REFERENCE_CONTINUATION\n' > "$dest/status.txt"
sha256sum -c "$dest/runner.sha256" > "$dest/prelaunch_runner_integrity.txt"
sha256sum -c "$dest/source.sha256" > "$dest/prelaunch_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/prelaunch_input_integrity.txt"
sha256sum -c "$dest/binary.sha256" > "$dest/prelaunch_binary_integrity.txt"
python3 - "$dest" <<'CHECK'
from pathlib import Path
import sys,shutil
text=Path('/proc/meminfo').read_text()
memory=dict((x.split(':',1)[0],x.split(':',1)[1].strip()) for x in text.splitlines())
assert int(memory['MemAvailable'].split()[0])*1024>128*2**30
assert shutil.disk_usage(sys.argv[1]).free>8*2**30
(Path(sys.argv[1])/'host_memory_at_launch.txt').write_text(text)
CHECK
# Slurm owns one typed GRES. Do not reject jobs on unallocated physical GPUs
# merely because nvidia-smi can enumerate the node's other devices.
nvidia-smi --query-gpu=index,uuid,name,memory.total,memory.used --format=csv > "$dest/device.txt"
printf 'job=%s\nnode=%s\nCUDA_VISIBLE_DEVICES=%s\n' "$SLURM_JOB_ID" "$(hostname)" "${CUDA_VISIBLE_DEVICES:-unset}" > "$dest/allocation.txt"
# Both H200 and any later lifetime-split candidate must acquire this exact
# atomic claim. Preserve it after termination for accepted-state audit.
if ! mkdir "$claim"; then
    printf 'REFUSED_EXISTING_SHARED_STEP6_EXECUTION_CLAIM\n' > "$dest/status.txt"
    exit 4
fi
claim_owned=1
cp "$dest/execution_contract.json" "$claim/execution_contract.json"
printf 'job=%s\nnode=%s\noutput=%s\n' "$SLURM_JOB_ID" "$(hostname)" "$dest/run" > "$claim/owner.txt"
mkdir "$dest/run"
printf 'RUNNING_H200_REFERENCE_CACHE_AND_TWO_CONTINUOUS_NATIVE_STEPS\n' > "$dest/status.txt"
set +e
timeout --signal=TERM --kill-after=30s 1800s mpirun -np 4 "$donor/continuation" \
    "$root/material_full_initial_v2_h200_mpi4/run" \
    "$root/material_full_positive_rule_v1_ada_mpi4" \
    "$root/material_full_positive_transition_v2_ada_mpi4/immutable_payload.sha256" \
    "$root/material_full_positive_transition_v2_ada_mpi4/run/continuous_step5/accepted_r" \
    "$root/material_cell_condensation_structure_v1_mpi4/run" \
    "$root/material_cell_condensation_factor_v1_mpi4/run" \
    "$root/material_cell_condensation_actual_rhs_v1_ada_mpi4/run" \
    "$donor/reference_source_binding.txt" "$dest/run" > "$dest/native.log" 2>&1
native_result=$?
set -e
printf '%s\n' "$native_result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/final_input_integrity.txt"
sha256sum -c "$dest/binary.sha256" > "$dest/final_binary_integrity.txt"
sha256sum -c "$dest/runner.sha256" > "$dest/final_runner_integrity.txt"
if test "$native_result" = 0; then
    printf 'TWO_NATIVE_STEPS_COMPLETE_LONG_RELEASE_NOT_AUTHORIZED_BY_GATES\n' > "$dest/status.txt"
else
    printf 'H200_REFERENCE_CACHE_CONTINUATION_REJECTED_%s\n' "$native_result" > "$dest/status.txt"
fi
exit "$native_result"
