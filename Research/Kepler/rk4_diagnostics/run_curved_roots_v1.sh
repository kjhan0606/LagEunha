#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
dest=$root/rk4_diagnostics/curved_roots_$version
case "$(hostname)" in LagEunha|lageunha) ;; *) printf 'Use authorized lageunha\n' >&2;exit 2;; esac
mkdir "$dest"
source "$root/run_env.sh"
export OMP_NUM_THREADS=2 I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
for h in gfs_material_map12.h gfs_material_p3.h gfs_material_p3_radial_input.h gfs_material_triangle.h gfs_moment.h; do cp "$root/rk4_diagnostics/material_p3_finite_v2/$h" "$dest/"; done
for f in material_curved_packet.h material_curved_bounds.h material_curved_compiled.h material_curved_roots.h test_curved_roots_v1.cpp; do cp "$root/rk4_diagnostics/$f" "$dest/"; done
packet=$root/rk4_diagnostics/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
sites=$root/rk4_diagnostics/curved_map_algebra_v2/advected_sites.csv
sha256sum "$packet" "$sites" > "$dest/inputs.sha256"
mpicxx -O2 -std=c++17 -fopenmp -fno-fast-math -ffp-contract=off -I"$dest" "$dest/test_curved_roots_v1.cpp" -o "$dest/contract.exe" > "$dest/build.log" 2>&1
for ranks in 1 4; do
    mkdir "$dest/mpi$ranks"
    mpirun -np "$ranks" "$dest/contract.exe" "$packet" "$sites" "$dest/mpi$ranks" > "$dest/mpi$ranks/contract.json" 2> "$dest/mpi$ranks/stderr.log"
done
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
cmp "$dest/mpi1/angular_partition.csv" "$dest/mpi4/angular_partition.csv"
printf 'AUDIT_COMPLETE\n' > "$dest/status.txt"
