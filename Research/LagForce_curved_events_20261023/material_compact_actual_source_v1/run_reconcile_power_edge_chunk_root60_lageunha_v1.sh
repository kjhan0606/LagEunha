#!/bin/bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
dest=$root/material_compact_actual_source_v1
transition=$root/material_full_positive_transition_v2_ada_mpi4
test "$#" = 2
first=$1
last=$2
[[ "$first" =~ ^[0-9]+$ && "$last" =~ ^[0-9]+$ ]]
(( last > first && last - first <= 50000 ))
chunk=$dest/power_edge_chunk_${first}_${last}_h200sm90_root60
if ! test -f "$chunk/native_exit_code.txt" ||
   test "$(cat "$chunk/native_exit_code.txt")" != 0; then
  chunk=$dest/power_edge_chunk_${first}_${last}_lageunha_root60
fi
presence=$dest/source_edge_presence_v3
run=$dest/power_edge_reconcile_${first}_${last}_root60_v1
test "$(hostname)" = LagEunha
test "$(cat "$chunk/native_exit_code.txt")" = 0
test "$(cat "$presence/native_exit_code.txt")" = 0
source "$transition/run_env.sh"
export I_MPI_CXX=g++ OMP_NUM_THREADS=1 OMP_PROC_BIND=false
cd "$dest"
sha256sum -c "$chunk/source_edge_wires.sha256" > /dev/null
sha256sum -c "$presence/presence.sha256" > /dev/null
mkdir "$run"
sha256sum reconcile_power_edge_chunk.cpp run_reconcile_power_edge_chunk_root60_lageunha_v1.sh \
  "$root/material_compact_rezone_parallel.cpp" \
  "$root/material_compact_overlap_events.h" \
  "$transition/Exam/libexam.a" \
  "$chunk"/source_edge_wires_rank{0,1,2,3}.csv \
  "$chunk"/source_edge_manifest_rank{0,1,2,3}.csv \
  "$presence/external_edge.bin" > "$run/source.sha256"
mpicxx -std=gnu++17 -fext-numeric-literals -O2 -ffunction-sections -fdata-sections \
  -fopenmp -fno-fast-math -ffp-contract=off -DXYZDBL -DGOTPM -DUSE_CUDA \
  -I /gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src/Exam \
  -I /opt/ohpc/pub/cuda/13.0.2/include reconcile_power_edge_chunk.cpp \
  -Wl,--gc-sections -L "$transition/Exam" -lexam \
  -L /opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lquadmath \
  -o "$run/reconcile_power_edge_chunk"
sha256sum "$run/reconcile_power_edge_chunk" > "$run/binary.sha256"
set +e
timeout --signal=TERM --kill-after=30s 600s \
  "$run/reconcile_power_edge_chunk" "$chunk" \
  "$presence/external_edge.bin" "$run" 642787 \
  </dev/null > "$run/native.log" 2>&1
result=$?
set -e
printf '%s\n' "$result" > "$run/native_exit_code.txt"
sha256sum -c "$run/source.sha256" > "$run/final_source_check.txt"
if test "$result" = 0; then
  sha256sum "$run"/canonical_owner{0,1,2,3}.csv > "$run/canonical.sha256"
fi
tail -8 "$run/native.log"
exit "$result"
