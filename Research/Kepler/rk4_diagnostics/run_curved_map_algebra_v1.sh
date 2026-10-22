#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
packet=${2:-$root/rk4_diagnostics/material_p3_finite_v1/mpi1_cpu/curved_map_control.bin}
dest=$root/rk4_diagnostics/curved_map_algebra_$version
case "$(hostname)" in LagEunha|lageunha) ;; *) printf 'Use authorized lageunha\n' >&2;exit 2;; esac
mkdir "$dest"
source "$root/run_env.sh"
export OMP_NUM_THREADS=1
headers=(gfs_material_map12.h gfs_material_p3.h gfs_material_p3_radial_input.h gfs_material_triangle.h gfs_moment.h gfs_moment_density.h gfs_moment_initial.h gfs_dg_euler.h gfs_enriched_euler.h gfs_ale_flux.h)
for h in "${headers[@]}"; do cp "$root/rk4_diagnostics/material_p3_finite_v1/$h" "$dest/"; done
cp "$root/rk4_diagnostics/material_curved_packet.h" "$root/rk4_diagnostics/test_curved_map_algebra_v1.cpp" "$root/rk4_diagnostics/test_enriched_euler_v1.cpp" "$dest/"
sha256sum "$packet" > "$dest/packet.sha256"
g++ -O2 -std=c++17 -I"$dest" "$dest/test_curved_map_algebra_v1.cpp" -o "$dest/contract.exe" > "$dest/build.log" 2>&1
"$dest/contract.exe" "$packet" "$root/rk4_diagnostics/dg_projection_v2_p2/fine_fresh" "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
printf 'AUDIT_COMPLETE\n' > "$dest/status.txt"
