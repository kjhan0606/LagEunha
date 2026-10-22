#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
ranks=${2:-4}
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
case "$ranks" in 1|4) ;; *) exit 2;; esac
dest=$root/rk4_diagnostics/shared_clipped_target_${version}_mpi${ranks}
mkdir "$dest"
source "${GFS_CURVED_FROZEN_ENV:-$root/run_env.sh}"
export OMP_NUM_THREADS=4 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
retained=$root/rk4_diagnostics/shared_target_v3_mpi4
bounds=$root/rk4_diagnostics/shared_real_bound_v1_mpi4
pressure=$root/rk4_diagnostics/geometric_pressure_v2_mpi4
coverage=$root/rk4_diagnostics/shared_boundary_witness_v5_mpi1
wholegroup=$root/rk4_diagnostics/shared_atlas_v4_mpi4
for f in "$retained/"*.h; do cp "$f" "$dest/"; done
cp "$retained/material_shared_duffy.h" "$dest/retained_duffy.h"
for f in test_shared_trace_factory_v1.cpp test_geometric_pressure_v1.cpp test_shared_function_v1.cpp material_shared_real_bound.h material_shared_duffy.h material_shared_target_geometry.h material_shared_target_integral.h material_shared_winding.h material_shared_retained.h material_shared_publication.h material_shared_pullback.h material_shared_atlas.h material_shared_clip.h material_shared_reference80_ball.h material_shared_reference80_integral.h test_shared_clipped_target_v1.cpp; do cp "$root/rk4_diagnostics/$f" "$dest/"; done
sha256sum -c "$retained/source.sha256" > "$dest/retained_source_proof_integrity.txt"
sha256sum -c "$retained/output.sha256" > "$dest/retained_output_proof_integrity.txt"
sha256sum -c "$bounds/output.sha256" > "$dest/retained_real_proof_integrity.txt"
sha256sum -c "$coverage/output.sha256" > "$dest/retained_coverage_integrity.txt"
sha256sum -c "$wholegroup/output.sha256" > "$dest/retained_whole_group_integrity.txt"
packet=$root/rk4_diagnostics/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
saved=$root/rk4_diagnostics/shared_trace_rhs_v1_mpi4
sites=$root/rk4_diagnostics/curved_map_algebra_v2/advected_sites.csv
weights=$root/rk4_diagnostics/curved_events_v9/site_weights.csv
previous=$root/rk4_diagnostics/shared_clipped_target_v7_mpi4
tau=$root/rk4_diagnostics/curved_whole_target_v3_mpi4/whole_target_ledger.csv
sha256sum "$dest/"*.h "$dest/"*.cpp > "$dest/source.sha256"
sha256sum "$packet" "$saved/shared_state.bin" "$sites" "$weights" "$retained/source_target_contributions.csv" "$retained/frozen_geometric_EOS_context.csv" "$retained/actual_target_"*.csv "$bounds/source_real_bounds.csv" "$coverage/"*.csv "$wholegroup/"*.csv "$tau" "$previous/"*.csv "$pressure/geometric_cuda.o" > "$dest/retained_input.sha256"
printf '%s\n' "${GFS_CURVED_LAUNCH_ENVELOPE:-missing}" > "$dest/runner_envelope.txt"
set +e
mpicxx -O2 -std=c++17 -fext-numeric-literals -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_shared_clipped_target_v1.cpp" "$pressure/geometric_cuda.o" -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lquadmath -o "$dest/contract.exe" > "$dest/build.log" 2>&1
build=$?
set -e
printf '%s\n' "$build" > "$dest/build_exit_code.txt"
if test "$build" != 0; then printf 'BUILD_FAILED\n' > "$dest/status.txt"; exit 1;fi
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
set +e
mpirun -np "$ranks" "$dest/contract.exe" "$packet" "$saved" "$sites" "$weights" "$retained" "$bounds" "$coverage" "$wholegroup" "$tau" "$previous" "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/retained_input.sha256" > "$dest/final_input_integrity.txt"
shopt -s nullglob
outputs=("$dest/"*.json "$dest/"*.csv "$dest/"*.bin)
if test "${#outputs[@]}" -gt 0; then sha256sum "${outputs[@]}" > "$dest/output.sha256"; fi
if test "$result" = 0; then printf 'WHOLE_ORIGINAL_TARGET_ALL12_PUBLISHED\n' > "$dest/status.txt"; else printf 'ACTUAL_CLIPPED_WHOLE_TARGET_REJECTED_EXIT_%s\n' "$result" > "$dest/status.txt"; fi
exit "$result"
