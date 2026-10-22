#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
ranks=${2:-4}
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
case "$ranks" in 1|4) ;; *) exit 2;; esac
dest=$root/rk4_diagnostics/shared_target_${version}_mpi${ranks}
mkdir "$dest"
mkdir "$dest/duffy_preflight"
source "${GFS_CURVED_FROZEN_ENV:-$root/run_env.sh}"
export OMP_NUM_THREADS=4 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
proof=$root/rk4_diagnostics/shared_function_v2_mpi4
pressure=$root/rk4_diagnostics/geometric_pressure_v2_mpi4
for f in "$proof/"*.h; do cp "$f" "$dest/"; done
for f in test_shared_trace_factory_v1.cpp test_geometric_pressure_v1.cpp test_shared_function_v1.cpp material_shared_real_bound.h material_shared_duffy.h material_shared_target_geometry.h material_shared_target_integral.h material_shared_winding.h material_shared_retained.h material_shared_publication.h test_shared_duffy_v1.cpp test_shared_target_v1.cpp; do cp "$root/rk4_diagnostics/$f" "$dest/"; done
for f in material_shared_function.h material_shared_native_ball.h material_forward_ball.h material_shared_trace_batch.h; do cmp "$dest/$f" "$proof/$f"; done
sha256sum -c "$proof/source.sha256" > "$dest/retained_source_proof_integrity.txt"
sha256sum -c "$proof/retained_input.sha256" > "$dest/retained_input_proof_integrity.txt"
sha256sum -c "$proof/output.sha256" > "$dest/retained_output_proof_integrity.txt"
packet=$root/rk4_diagnostics/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
saved=$root/rk4_diagnostics/shared_trace_rhs_v1_mpi4
sites=$root/rk4_diagnostics/curved_map_algebra_v2/advected_sites.csv
weights=$root/rk4_diagnostics/curved_events_v9/site_weights.csv
normalizers=$root/rk4_diagnostics/curved_whole_target_v3_mpi4/whole_target_ledger.csv
sha256sum "$dest/"*.h "$dest/"*.cpp > "$dest/source.sha256"
sha256sum "$packet" "$saved/shared_state.bin" "$sites" "$weights" "$normalizers" "$pressure/finite/"frozen_physical_rank*.bin "$pressure/geometric_cuda.o" > "$dest/retained_input.sha256"
printf '%s\n' "${GFS_CURVED_LAUNCH_ENVELOPE:-missing}" > "$dest/runner_envelope.txt"
for stem in shared_duffy shared_target; do
  set +e
  mpicxx -O2 -std=c++17 -fext-numeric-literals -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_${stem}_v1.cpp" "$pressure/geometric_cuda.o" -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lquadmath -o "$dest/${stem}.exe" > "$dest/${stem}_build.log" 2>&1
  build=$?
  set -e
  printf '%s\n' "$build" > "$dest/${stem}_build_exit_code.txt"
  if test "$build" != 0; then printf 'BUILD_FAILED_%s\n' "$stem" > "$dest/status.txt"; exit 1;fi
done
sha256sum "$dest/"*.exe > "$dest/binary.sha256"
set +e
mpirun -np "$ranks" "$dest/shared_duffy.exe" "$packet" "$saved" "$sites" "$weights" "$dest/duffy_preflight" > "$dest/duffy_preflight/contract.json" 2> "$dest/duffy_preflight/stderr.log"
preflight=$?
set -e
printf '%s\n' "$preflight" > "$dest/duffy_preflight/native_exit_code.txt"
if test "$preflight" != 0; then printf 'DUFFY_PREFLIGHT_REJECTED\n' > "$dest/status.txt"; exit "$preflight";fi
set +e
mpirun -np "$ranks" "$dest/shared_target.exe" "$packet" "$saved" "$pressure" "$sites" "$weights" "$normalizers" "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/retained_input.sha256" > "$dest/final_input_integrity.txt"
sha256sum "$dest/"*.json "$dest/"*.csv "$dest/"*.bin "$dest/duffy_preflight/"*.json "$dest/duffy_preflight/"*.csv "$dest/duffy_preflight/"*.bin > "$dest/output.sha256"
if test "$result" = 0; then printf 'WHOLE_TARGET_ALL12_PUBLISHED\n' > "$dest/status.txt"; else printf 'WHOLE_TARGET_TRANSACTION_REJECTED_EXIT_%s\n' "$result" > "$dest/status.txt"; fi
exit "$result"
