#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
src=/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src
version=${1:-v1}
dest=$root/rk4_diagnostics/material_p3_hybrid_$version
case "$(hostname)" in LagEunha|lageunha) ;; *) printf 'Use authorized lageunha\n' >&2;exit 2;; esac
mkdir "$dest"
source "$root/run_env.sh"
export OMP_NUM_THREADS=1 I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
cd "$src"
headers=(gfs_material_p3.h gfs_material_p3_radial_input.h gfs_material_map.h gfs_material_triangle.h gfs_material_trace.h gfs_moment.h gfs_moment_density.h gfs_moment_initial.h gfs_dg_euler.h gfs_enriched_euler.h gfs_ale_flux.h gfs_ale_fan.h)
for h in "${headers[@]}"; do sha256sum "Exam/$h" >> "$dest/source.sha256"; cp "Exam/$h" "$dest/"; done
for f in test_enriched_euler_v1.cpp test_material_shared_v1.cpp test_material_hybrid_v1.cpp material_p3_mesh.h material_p3_trace.h test_material_p3_hybrid_v1.cpp; do cp "$root/rk4_diagnostics/$f" "$dest/"; done
mpicxx -O2 -std=c++17 -I"$dest" "$dest/test_material_p3_hybrid_v1.cpp" -o "$dest/contract.exe" > "$dest/build.log" 2>&1
for ranks in 1 4; do
    mkdir "$dest/mpi$ranks"
    mpirun -np "$ranks" "$dest/contract.exe" "$root/rk4_diagnostics/dg_projection_v2_p2/fine_fresh" "$dest/mpi$ranks" ${2:+"$2"} > "$dest/mpi$ranks/contract.json" 2> "$dest/mpi$ranks/stderr.log"
done
sha256sum -c "$dest/source.sha256" > "$dest/source_unchanged.txt"
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
printf 'AUDIT_COMPLETE\n' > "$dest/status.txt"
