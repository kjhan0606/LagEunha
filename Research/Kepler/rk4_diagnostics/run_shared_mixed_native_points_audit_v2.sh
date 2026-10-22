#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
diag=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
cd "$diag"
run=$diag/shared_mixed_native_v2_mpi1
dest=$diag/shared_mixed_native_points_audit_v2
test "$(cat "$run/status.txt")" = MIXED_NATIVE_READONLY_EXIT_0
test "$(cat "$run/native_exit_code.txt")" = 0
test ! -e "$dest/native_exit_code.txt"
test ! -e "$dest/contract.json"
sha256sum -c "$dest/source.sha256" > "$dest/source_integrity.txt"
sha256sum -c "$dest/binary.sha256" > "$dest/binary_integrity.txt"
sha256sum -c "$run/output.sha256" > "$dest/native_output_integrity.txt"
sha256sum "$run/source.sha256" "$run/input.sha256" "$run/output.sha256" "$run/native_source_packets_rank0.bin" > "$dest/input.sha256"
source "$GFS_CURVED_FROZEN_ENV"
export OMP_NUM_THREADS=16 OMP_PROC_BIND=false
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY
printf '%s\n' "$GFS_CURVED_LAUNCH_ENVELOPE" > "$dest/runner_envelope.txt"
set +e
"$dest/audit.exe" "$run" "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/input.sha256" > "$dest/final_input_integrity.txt"
sha256sum "$dest/contract.json" "$dest/exact_native_prefix_accounting.csv" > "$dest/output.sha256"
printf 'NATIVE_INTEGER_AUDIT_EXIT_%s\n' "$result" > "$dest/status.txt"
exit "$result"
