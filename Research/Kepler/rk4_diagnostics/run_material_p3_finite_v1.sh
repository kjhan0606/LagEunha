#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
src=/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src
version=${1:-v1}
dest=$root/rk4_diagnostics/material_p3_finite_$version
case "$(hostname)" in LagEunha|lageunha) ;; *) printf 'Use authorized lageunha\n' >&2;exit 2;; esac
active=$(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader)
if test -n "$active"; then printf 'GPU compute process present; no launch\n' >&2;exit 3;fi
mkdir "$dest"
source "$root/run_env.sh"
export OMP_NUM_THREADS=2 I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++ CUDA_VISIBLE_DEVICES=0
unset I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
cd "$src"
headers=(gfs_material_map12.h gfs_material_p3.h gfs_material_p3_radial_input.h gfs_material_map.h gfs_material_triangle.h gfs_material_trace.h gfs_moment.h gfs_moment_density.h gfs_moment_initial.h gfs_dg_euler.h gfs_enriched_euler.h gfs_ale_flux.h gfs_ale_fan.h)
for h in "${headers[@]}"; do sha256sum "Exam/$h" >> "$dest/source.sha256"; cp "Exam/$h" "$dest/"; done
files=(test_enriched_euler_v1.cpp test_material_shared_v1.cpp test_material_hybrid_v1.cpp material_p3_mesh.h material_p3_trace.h test_material_p3_hybrid_v1.cpp material_p3_finite_batch.h material_p3_finite_cuda.cu test_material_p3_finite_v1.cpp)
for f in "${files[@]}"; do cp "$root/rk4_diagnostics/$f" "$dest/"; done
nvcc -O2 -std=c++17 -I"$dest" -c "$dest/material_p3_finite_cuda.cu" -o "$dest/finite_cuda.o" > "$dest/build_cuda.log" 2>&1
mpicxx -O2 -std=c++17 -fopenmp -I"$dest" "$dest/test_material_p3_finite_v1.cpp" "$dest/finite_cuda.o" -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -o "$dest/contract.exe" > "$dest/build.log" 2>&1
nvidia-smi --query-gpu=index,name,memory.total,memory.used --format=csv > "$dest/device.txt"
seed=$root/rk4_diagnostics/material_p3_hybrid_v4/mpi1/level1_Kepler_state.bin
sha256sum "$seed" > "$dest/seed.sha256"
for backend in cpu cuda; do
    for ranks in 1 4; do
        mkdir "$dest/mpi${ranks}_$backend"
        mpirun -np "$ranks" "$dest/contract.exe" "$root/rk4_diagnostics/dg_projection_v2_p2/fine_fresh" "$seed" "$dest/mpi${ranks}_$backend" "$backend" ${2:+"$2"} > "$dest/mpi${ranks}_$backend/contract.json" 2> "$dest/mpi${ranks}_$backend/stderr.log"
    done
done
sha256sum -c "$dest/source.sha256" > "$dest/source_unchanged.txt"
sha256sum "$dest/contract.exe" > "$dest/binary.sha256"
printf 'AUDIT_COMPLETE\n' > "$dest/status.txt"
