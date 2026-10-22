#!/usr/bin/env bash
#SBATCH --job-name=shared-evolve
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=00:20:00
set -euo pipefail
dest=$1
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
source "$dest/run_env.sh"
export OMP_NUM_THREADS=4 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
export GFS_EVOLUTION_RETAINED_MANIFEST=$dest/retained_input.sha256
packet=$root/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
saved=$root/shared_trace_rhs_v1_mpi4/shared_state.bin
printf 'BUILDING\n' > "$dest/status.txt"
sha256sum -c "$dest/source.sha256" > "$dest/prelaunch_source_integrity.txt"
sha256sum -c "$dest/retained_input.sha256" > "$dest/prelaunch_input_integrity.txt"
nvidia-smi --query-gpu=index,name,memory.total,memory.used --format=csv > "$dest/device.txt"
set +e
nvcc -O2 -std=c++17 --fmad=false -I"$dest" -c "$dest/material_geometric_pressure_persistent.cu" -o "$dest/geometric_cuda.o" > "$dest/build_cuda.log" 2>&1
cuda_build=$?
if test "$cuda_build" = 0; then
    mpicxx -O2 -std=c++17 -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_shared_evolution_v1.cpp" "$dest/geometric_cuda.o" -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -o "$dest/contract.exe" > "$dest/build.log" 2>&1
    host_build=$?
else host_build=1; fi
set -e
printf '%s\n' "$cuda_build" > "$dest/build_cuda_exit_code.txt"
printf '%s\n' "$host_build" > "$dest/build_exit_code.txt"
if test "$cuda_build" != 0 || test "$host_build" != 0; then printf 'BUILD_FAILED\n' > "$dest/status.txt"; exit 1; fi
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
printf 'RUNNING\n' > "$dest/status.txt"
set +e
mpirun -np 4 "$dest/contract.exe" "$packet" "$saved" "$root" "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/retained_input.sha256" > "$dest/final_input_integrity.txt"
if test "$result" = 0; then
    sha256sum "$dest/"*.bin "$dest/contract.json" "$dest/accepted_first/step.json" "$dest/uninterrupted_second/step.json" "$dest/restarted_second/step.json" > "$dest/output.sha256"
    printf 'ACTUAL_RK4_RESTART_PASS\n' > "$dest/status.txt"
else printf 'ACTUAL_RK4_REJECTED_EXIT_%s\n' "$result" > "$dest/status.txt"; fi
exit "$result"
