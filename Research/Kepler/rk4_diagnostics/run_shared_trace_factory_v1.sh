#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
ranks=${2:-4}
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
case "$ranks" in 1|4) ;; *) exit 2;; esac
dest=$root/rk4_diagnostics/shared_trace_factory_${version}_mpi${ranks}
mkdir "$dest"
active=$(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader)
if test -n "$active"; then printf '%s\n' "$active" > "$dest/gpu_busy.txt"; printf 'GPU_BUSY_NO_LAUNCH\n' > "$dest/status.txt"; exit 3;fi
source "${GFS_CURVED_FROZEN_ENV:-$root/run_env.sh}"
export OMP_NUM_THREADS=4 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++ CUDA_VISIBLE_DEVICES=0
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
saved=$root/rk4_diagnostics/material_p3_finite_v2
for f in "$saved/"gfs*.h; do cp "$f" "$dest/"; done
for f in material_curved_packet.h material_shared_trace_batch.h material_shared_trace_factory.h material_shared_trace_cuda.cu test_shared_trace_factory_v1.cpp; do cp "$root/rk4_diagnostics/$f" "$dest/"; done
packet=$root/rk4_diagnostics/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
sha256sum "$dest/"*.h "$dest/"*.cpp "$dest/"*.cu > "$dest/source.sha256"
sha256sum "$packet" > "$dest/retained_input.sha256"
printf '%s\n' "${GFS_CURVED_LAUNCH_ENVELOPE:-missing}" > "$dest/runner_envelope.txt"
nvidia-smi --query-gpu=index,name,memory.total,memory.used --format=csv > "$dest/device.txt"
set +e
nvcc -O2 -std=c++17 --fmad=false -I"$dest" -c "$dest/material_shared_trace_cuda.cu" -o "$dest/shared_cuda.o" > "$dest/build_cuda.log" 2>&1
cuda_build=$?
mpicxx -O2 -std=c++17 -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_shared_trace_factory_v1.cpp" "$dest/shared_cuda.o" -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -o "$dest/contract.exe" > "$dest/build.log" 2>&1
host_build=$?
set -e
printf '%s\n' "$cuda_build" > "$dest/build_cuda_exit_code.txt"
printf '%s\n' "$host_build" > "$dest/build_exit_code.txt"
if test "$cuda_build" != 0 || test "$host_build" != 0; then printf 'BUILD_FAILED\n' > "$dest/status.txt"; exit 1;fi
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
set +e
mpirun -np "$ranks" "$dest/contract.exe" "$packet" "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/retained_input.sha256" > "$dest/final_input_integrity.txt"
sha256sum "$dest/"*.json "$dest/"*.csv "$dest/"*.bin > "$dest/output.sha256"
if test "$result" = 0; then printf 'FACTORY_INITIAL_ONLY_PASS_RHS_FINITE_CONTROL_PENDING\n' > "$dest/status.txt"; else printf 'FACTORY_REJECTED_EXIT_%s\n' "$result" > "$dest/status.txt"; fi
exit "$result"
