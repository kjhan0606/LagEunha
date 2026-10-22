#!/usr/bin/env bash
#SBATCH --job-name=full-material-rhs
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=128G
#SBATCH --time=00:30:00
set -euo pipefail
dest=$1
source "$dest/run_env.sh"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
sha256sum -c "$dest/source.sha256" > "$dest/source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/input_integrity.txt"
sha256sum -c "$dest/prebuild_binary.sha256" > "$dest/binary_integrity.txt"
nvidia-smi --query-gpu=index,name,uuid,memory.total,memory.used --format=csv > "$dest/devices.txt"
printf 'EVALUATING_ACTUAL_FULL_PHYSICAL_RHS\n' > "$dest/status.txt"
set +e
mpirun -np 4 "$dest/prebuild.exe" /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/shared_full_initialization_v3_h200_mpi4 "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/final_input_integrity.txt"
sha256sum "$dest/"*.bin "$dest/"*.json "$dest/"*.csv > "$dest/output.sha256"
if test "$result" = 0; then printf 'ACTUAL_FULL_RHS_PASS_NOT_YET_TIMESTEPPED\n' > "$dest/status.txt"; else printf 'FULL_PHYSICAL_RHS_REJECTED_EXIT_%s\n' "$result" > "$dest/status.txt"; fi
exit "$result"
