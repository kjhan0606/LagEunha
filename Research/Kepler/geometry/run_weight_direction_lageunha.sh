#!/bin/bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
dest=$root/material_compact_actual_source_v1
input=$dest/power_v25
run=$dest/weight_direction_v1
test "$(hostname)" = LagEunha
test "$(cat "$input/native_exit_code.txt")" = 0
cd "$dest"
sha256sum -c "$input/gauss_area_residual.sha256" > weight_area_input_check.txt
sha256sum -c "$input/power_area_jacobian_edges.sha256" > weight_edge_input_check.txt
mkdir "$run"
sha256sum solve_weight_direction.cpp "$input/gauss_area_residual.csv" \
  "$input/power_area_jacobian_edges.csv" > "$run/source_and_inputs.sha256"
g++ -std=gnu++17 -O2 -fno-fast-math -ffp-contract=off \
  solve_weight_direction.cpp -o "$run/solve_weight_direction"
sha256sum "$run/solve_weight_direction" > "$run/binary.sha256"
TIMEFORMAT='elapsed_seconds=%3R'
set +e
time timeout --signal=TERM --kill-after=10s 300s \
  "$run/solve_weight_direction" "$input/gauss_area_residual.csv" \
  "$input/power_area_jacobian_edges.csv" "$run/weights.csv" \
  > "$run/native.log" 2>&1
result=$?
set -e
printf '%s\n' "$result" > "$run/native_exit_code.txt"
sha256sum -c "$run/source_and_inputs.sha256" > "$run/final_source_check.txt"
if test "$result" = 0; then sha256sum "$run/weights.csv" > "$run/weights.sha256"; fi
cat "$run/native.log"
exit "$result"
