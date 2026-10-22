#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
src=/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src
version=${1:-v1}
dest=$root/rk4_diagnostics/material_remap_$version
case "$(hostname)" in LagEunha|lageunha) ;; *) printf 'Use authorized lageunha\n' >&2;exit 2;; esac
mkdir "$dest"
source "$root/run_env.sh"
export OMP_NUM_THREADS=1
cd "$src"
sha256sum Exam/gfs_moment.h Exam/gfs_moment_density.h Exam/gfs_moment_initial.h Exam/gfs_dg_euler.h Exam/gfs_enriched_euler.h Exam/gfs_ale_flux.h > "$dest/source.sha256"
cp Exam/gfs_moment.h Exam/gfs_moment_density.h Exam/gfs_moment_initial.h Exam/gfs_dg_euler.h Exam/gfs_enriched_euler.h Exam/gfs_ale_flux.h "$dest/"
cp "$root/rk4_diagnostics/test_enriched_euler_v1.cpp" "$root/rk4_diagnostics/test_material_remap_v1.cpp" "$dest/"
g++ -O2 -std=c++17 -I"$dest" "$dest/test_material_remap_v1.cpp" -o "$dest/contract.exe" > "$dest/build.log" 2>&1
"$dest/contract.exe" "$root/rk4_diagnostics/dg_projection_v2_p2/fine_fresh" > "$dest/contract.json" 2> "$dest/stderr.log"
sha256sum -c "$dest/source.sha256" > "$dest/source_unchanged.txt"
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
printf 'AUDIT_COMPLETE\n' > "$dest/status.txt"
