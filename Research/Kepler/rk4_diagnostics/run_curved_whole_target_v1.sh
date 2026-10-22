#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
ranks=${2:-4}
replay=${3:-}
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
case "$ranks" in 1|4) ;; *) exit 2;; esac
dest=$root/rk4_diagnostics/curved_whole_target_${version}_mpi${ranks}
guard=$root/rk4_diagnostics/curved_events_v9
saved=$root/rk4_diagnostics/curved_target_source_glue_v2_225
mkdir "$dest"
source "${GFS_CURVED_FROZEN_ENV:-$root/run_env.sh}"
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
for f in "$saved/"*.h; do cp "$f" "$dest/"; done
for f in test_curved_transfer_v1.cpp test_curved_complement_v1.cpp test_curved_branch_v1.cpp test_curved_original_oracle_v1.cpp test_curved_target_integral_v1.cpp material_target_analytic.h material_target_gluing.h material_whole_target_geometry.h material_whole_target.h test_curved_whole_target_v1.cpp; do cp "$root/rk4_diagnostics/$f" "$dest/"; done
packet=$root/rk4_diagnostics/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
sites=$root/rk4_diagnostics/curved_map_algebra_v2/advected_sites.csv
lower=$root/rk4_diagnostics/curved_target_lower_v1_225/target_lower.csv
sha256sum "$dest/"*.h "$dest/"*.cpp > "$dest/source.sha256"
sha256sum -c "$guard/inputs.sha256" > "$dest/input_integrity.txt"
sha256sum "$packet" "$sites" "$guard/site_weights.csv" "$lower" > "$dest/retained_inputs.sha256"
printf '%s\n' "${GFS_CURVED_LAUNCH_ENVELOPE:-missing}" > "$dest/runner_envelope.txt"
args=()
if test -n "$replay"; then
    test "$replay" = "$root/rk4_diagnostics/curved_whole_target_v2_mpi4"
    sha256sum -c "$replay/source.sha256" > "$dest/reused_source_integrity.txt"
    sha256sum -c "$replay/retained_inputs.sha256" > "$dest/reused_input_integrity.txt"
    sha256sum -c "$replay/output.sha256" > "$dest/reused_output_integrity.txt"
    for f in material_whole_target.h material_target_analytic.h material_target_general.h material_target_domain.h material_target_inverse.h material_forward_ball.h material_original_oracle.h test_curved_target_integral_v1.cpp; do cmp "$dest/$f" "$replay/$f"; done
    args+=("$replay/source_contributions.csv")
    sha256sum "$replay/source_contributions.csv" "$replay/"positive_rank*.bin "$replay/"sectors_rank*.csv > "$dest/retained_source_proofs.sha256"
fi
set +e
mpicxx -O2 -std=c++17 -fext-numeric-literals -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_curved_whole_target_v1.cpp" -lquadmath -o "$dest/contract.exe" > "$dest/build.log" 2>&1
build_result=$?
set -e
printf '%s\n' "$build_result" > "$dest/build_exit_code.txt"
if test "$build_result" != 0; then printf 'BUILD_FAILED_EXIT_%s\n' "$build_result" > "$dest/status.txt"; exit "$build_result"; fi
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
set +e
mpirun -np "$ranks" "$dest/contract.exe" "$packet" "$sites" "$guard" "$lower" "$dest" "${args[@]}" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/retained_inputs.sha256" > "$dest/final_retained_input_integrity.txt"
sha256sum "$dest/"*.csv "$dest/"*.json "$dest/"*.bin > "$dest/output.sha256"
if test "$result" = 0; then printf 'ONE_TARGET_REFERENCE_PACKET_PUBLISHED_NOT_REMAP_OR_EVOLUTION\n' > "$dest/status.txt"; else printf 'WHOLE_TARGET_REJECTED_EXIT_%s\n' "$result" > "$dest/status.txt"; fi
exit "$result"
