#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
source_dir=/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src/Exam
audit=$root/enriched_contract_v1
mkdir "$audit"
cd "$source_dir"
sha256sum gfs_moment.h gfs_moment_density.h gfs_moment_initial.h gfs_dg_euler.h gfs_ale_flux.h gfs_enriched_euler.h > "$audit/source.sha256"
cp gfs_moment.h gfs_moment_density.h gfs_moment_initial.h gfs_dg_euler.h gfs_ale_flux.h gfs_enriched_euler.h "$audit/"
cp "$root/test_enriched_euler_v1.cpp" "$audit/"
g++ -O2 -std=c++14 -I"$audit" "$audit/test_enriched_euler_v1.cpp" -o "$audit/test_enriched_euler.exe" > "$audit/build.log" 2>&1
sha256sum "$audit/test_enriched_euler.exe" > "$audit/binary.sha256"
"$audit/test_enriched_euler.exe" "$root/dg_projection_v2_p2/fine_fresh" > "$audit/results.json" 2> "$audit/stderr.log"
sha256sum -c "$audit/source.sha256" > "$audit/source_unchanged.txt"
printf 'CONTRACT_COMPLETED\n' > "$audit/status.txt"
