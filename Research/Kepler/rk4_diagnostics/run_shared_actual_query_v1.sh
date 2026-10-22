#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
query_source=${2:-301}
quotient_trace=${3:-0}
stable_reference=${4:-0}
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
[[ "$query_source" =~ ^[0-9]+$ ]]
test "$query_source" -lt 384
case "$quotient_trace" in 0|1) ;; *) exit 2;; esac
case "$stable_reference" in 0|1) ;; *) exit 2;; esac
originating=$root/rk4_diagnostics/shared_clipped_target_v8_mpi4
# This is a dependent diagnostic, never a second concurrent target call.
test -f "$originating/status.txt"
test "$(cat "$originating/native_exit_code.txt")" = 2
dest=$root/rk4_diagnostics/shared_actual_query_${version}_source${query_source}_mpi1
mkdir "$dest"
source "${GFS_CURVED_FROZEN_ENV:-$root/run_env.sh}"
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
for input in "$originating/"*.h "$originating/"*.cpp; do cp "$input" "$dest/"; done
cp "$root/rk4_diagnostics/test_shared_actual_query_v1.cpp" "$dest/"
if test "$stable_reference" = 1; then
    cp "$root/rk4_diagnostics/material_shared_reference80_ball.h" "$dest/"
fi
diagnostic_flags=()
if test "$quotient_trace" = 1; then
    # Instrument the failed precondition only, in a private copy. Every
    # original arithmetic expression and the original rejection stay intact.
    python3 - "$dest/material_forward_ball.h" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1])
original=p.read_text()
guard='need(lower>0,"forward quotient denominator separated");'
assert original.count(guard)==1
trace='if(!(lower>0)){GFS_FORWARD_QUOTIENT_FAILURE(a.n,a.e,b.n,b.e);}'
p.write_text(original.replace(guard,trace+guard))
PY
    diagnostic_flags=(-g -rdynamic -fno-omit-frame-pointer -fno-pie -no-pie)
fi
sha256sum -c "$originating/source.sha256" > "$dest/originating_source_integrity.txt"
sha256sum -c "$originating/output.sha256" > "$dest/originating_output_integrity.txt"
sha256sum -c "$originating/retained_input.sha256" > "$dest/originating_input_integrity.txt"
packet=$root/rk4_diagnostics/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
saved=$root/rk4_diagnostics/shared_trace_rhs_v1_mpi4
sites=$root/rk4_diagnostics/curved_map_algebra_v2/advected_sites.csv
weights=$root/rk4_diagnostics/curved_events_v9/site_weights.csv
retained=$root/rk4_diagnostics/shared_target_v3_mpi4
bounds=$root/rk4_diagnostics/shared_real_bound_v1_mpi4
pressure=$root/rk4_diagnostics/geometric_pressure_v2_mpi4
owner=$((query_source%4))
sha256sum "$dest/"*.h "$dest/"*.cpp > "$dest/source.sha256"
sha256sum "$packet" "$saved/shared_state.bin" "$sites" "$weights" "$retained/source_target_contributions.csv" "$retained/frozen_geometric_EOS_context.csv" "$bounds/source_real_bounds.csv" "$originating/all384_source_completion.csv" "$originating/actual_reference_panel_packets_rank${owner}.bin" "$pressure/geometric_cuda.o" > "$dest/retained_input.sha256"
printf '%s\n' "${GFS_CURVED_LAUNCH_ENVELOPE:-missing}" > "$dest/runner_envelope.txt"
set +e
mpicxx -O2 -std=c++17 -fext-numeric-literals -fopenmp -fno-fast-math -ffp-contract=off "${diagnostic_flags[@]}" -I"$dest" "$dest/test_shared_actual_query_v1.cpp" "$pressure/geometric_cuda.o" -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lquadmath -o "$dest/contract.exe" > "$dest/build.log" 2>&1
build=$?
set -e
printf '%s\n' "$build" > "$dest/build_exit_code.txt"
if test "$build" != 0; then printf 'BUILD_FAILED\n' > "$dest/status.txt"; exit 1; fi
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
set +e
mpirun -np 1 "$dest/contract.exe" "$packet" "$saved" "$sites" "$weights" "$retained" "$bounds" "$originating" "$query_source" "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/retained_input.sha256" > "$dest/final_input_integrity.txt"
sha256sum "$dest/"*.json "$dest/"*.csv "$dest/"*.bin > "$dest/output.sha256"
if test "$result" = 0; then printf 'ONE_ACTUAL_QUERY_DIAGNOSIS_RECORDED\n' > "$dest/status.txt"; else printf 'ACTUAL_QUERY_DIAGNOSTIC_FAILED_%s\n' "$result" > "$dest/status.txt"; fi
exit "$result"
