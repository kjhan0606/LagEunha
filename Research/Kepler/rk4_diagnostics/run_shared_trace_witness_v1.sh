#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
dest=$root/rk4_diagnostics/shared_trace_witness_${version}
saved=$root/rk4_diagnostics/shared_trace_rhs_v1_mpi4
mkdir "$dest"
source "${GFS_CURVED_FROZEN_ENV:-$root/run_env.sh}"
export OMP_NUM_THREADS=4 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
for f in "$saved/"*.h; do cp "$f" "$dest/"; done
cp "$saved/test_shared_trace_factory_v1.cpp" "$dest/"
cp "$root/rk4_diagnostics/test_shared_trace_witness_v1.cpp" "$dest/"
cp "$saved/shared_cuda.o" "$dest/"
packet=$root/rk4_diagnostics/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
sha256sum -c "$saved/source.sha256" > "$dest/retained_source_integrity.txt"
sha256sum -c "$saved/output.sha256" > "$dest/retained_output_integrity.txt"
sha256sum "$dest/"*.h "$dest/"*.cpp "$dest/shared_cuda.o" > "$dest/source.sha256"
sha256sum "$packet" "$saved/initial_state.bin" "$saved/shared_acceleration.csv" > "$dest/retained_input.sha256"
printf '%s\n' "${GFS_CURVED_LAUNCH_ENVELOPE:-missing}" > "$dest/runner_envelope.txt"
set +e
mpicxx -O2 -std=c++17 -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_shared_trace_witness_v1.cpp" "$dest/shared_cuda.o" -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -o "$dest/contract.exe" > "$dest/build.log" 2>&1
build=$?
set -e
printf '%s\n' "$build" > "$dest/build_exit_code.txt"
if test "$build" != 0; then printf 'BUILD_FAILED\n' > "$dest/status.txt"; exit "$build";fi
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
set +e
mpirun -np 4 "$dest/contract.exe" "$packet" "$saved" "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/retained_input.sha256" > "$dest/final_input_integrity.txt"
sha256sum "$dest/"*.json "$dest/"*.csv > "$dest/output.sha256"
printf 'READ_ONLY_WITNESS_EXIT_%s\n' "$result" > "$dest/status.txt"
exit "$result"
