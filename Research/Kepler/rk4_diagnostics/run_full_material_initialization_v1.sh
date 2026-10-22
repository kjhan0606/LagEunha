#!/usr/bin/env bash
#SBATCH --job-name=full-material-init
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:4
#SBATCH --mem=128G
#SBATCH --time=01:00:00
set -euo pipefail
dest=$1
source "$dest/run_env.sh"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
sha256sum -c "$dest/source.sha256" > "$dest/source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/input_integrity.txt"
nvidia-smi --query-gpu=index,name,uuid,memory.total,memory.used --format=csv > "$dest/devices.txt"
printf 'BUILDING\n' > "$dest/status.txt"
set +e
nvcc -O2 -arch=sm_90 -std=c++17 --fmad=false -Xptxas=-v -I"$dest" -c "$dest/material_stream_mass.cu" -o "$dest/stream_mass.o" > "$dest/build_cuda.log" 2>&1
gpu_build=$?
if test "$gpu_build" = 0; then
    mpicxx -O2 -std=c++17 -fopenmp -fno-fast-math -ffp-contract=off -I/opt/ohpc/pub/cuda/13.0.2/include -I"$dest" "$dest/initialize_full_material.cpp" "$dest/stream_mass.o" -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -o "$dest/initialize.exe" > "$dest/build.log" 2>&1
    host_build=$?
else host_build=1; fi
set -e
printf '%s\n' "$gpu_build" > "$dest/build_cuda_exit_code.txt"
printf '%s\n' "$host_build" > "$dest/build_exit_code.txt"
if test "$gpu_build" != 0 || test "$host_build" != 0; then printf 'BUILD_FAILED\n' > "$dest/status.txt"; exit 1; fi
sha256sum "$dest/initialize.exe" > "$dest/binary.sha256"
printf 'INITIALIZING_FULL_REFERENCE\n' > "$dest/status.txt"
set +e
mpirun -np 4 "$dest/initialize.exe" /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/shared_full_reference_v1/full_reference.bin "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/final_input_integrity.txt"
if test "$result" = 0; then
    sha256sum "$dest/"*.bin "$dest/"*.json > "$dest/output.sha256"
    printf 'FULL_REFERENCE_MASS_INITIALIZED_NOT_EVOLVED\n' > "$dest/status.txt"
else printf 'FULL_INITIALIZATION_REJECTED_EXIT_%s\n' "$result" > "$dest/status.txt"; fi
exit "$result"
