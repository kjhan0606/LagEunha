#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
dest=$root/rk4_diagnostics/curved_bounds_$version
case "$(hostname)" in LagEunha|lageunha) ;; *) printf 'Use authorized lageunha\n' >&2;exit 2;; esac
mkdir "$dest"
source "$root/run_env.sh"
export OMP_NUM_THREADS=1
headers=(gfs_material_map12.h gfs_material_p3.h gfs_material_p3_radial_input.h gfs_material_triangle.h gfs_moment.h)
for h in "${headers[@]}"; do cp "$root/rk4_diagnostics/material_p3_finite_v2/$h" "$dest/"; done
for f in material_curved_packet.h material_curved_bounds.h test_curved_bounds_v1.cpp; do cp "$root/rk4_diagnostics/$f" "$dest/"; done
packet=$root/rk4_diagnostics/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
sha256sum "$packet" > "$dest/packet.sha256"
g++ -O2 -std=c++17 -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_curved_bounds_v1.cpp" -o "$dest/contract.exe" > "$dest/build.log" 2>&1
if test "$#" -ge 2; then
    parent=$2
    sha256sum -c "$parent/packet.sha256" > "$dest/parent_packet_unchanged.txt"
    sha256sum "$parent/interval_bounds.csv" > "$dest/resumed_ledger.sha256"
    printf '%s\n' "$parent" > "$dest/resumed_parent.txt"
    "$dest/contract.exe" "$packet" "$dest" "$parent/interval_bounds.csv" > "$dest/contract.json" 2> "$dest/stderr.log"
else
    "$dest/contract.exe" "$packet" "$dest" > "$dest/contract.json" 2> "$dest/stderr.log"
fi
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
printf 'AUDIT_COMPLETE\n' > "$dest/status.txt"
