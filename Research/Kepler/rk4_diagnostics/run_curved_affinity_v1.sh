#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
dest=$root/rk4_diagnostics/curved_affinity_v1
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
mkdir "$dest"
source "$root/run_env.sh"
export OMP_NUM_THREADS=2 I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
cp "$root/rk4_diagnostics/probe_curved_affinity.cpp" "$dest/"
mpicxx -O2 -std=c++17 -fopenmp "$dest/probe_curved_affinity.cpp" -o "$dest/probe.exe" > "$dest/build.log" 2>&1
mpirun -np 4 "$dest/probe.exe" > "$dest/inherited.txt" 2> "$dest/inherited.stderr"
export OMP_PROC_BIND=false
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY
mpirun -np 4 "$dest/probe.exe" > "$dest/unbound.txt" 2> "$dest/unbound.stderr"
mkdir "$dest/transfer_mpi4"
packet=$root/rk4_diagnostics/material_p3_finite_v2/mpi1_cpu/curved_map_control.bin
sites=$root/rk4_diagnostics/curved_map_algebra_v2/advected_sites.csv
exe=$root/rk4_diagnostics/curved_transfer_v3/contract.exe
sha256sum "$exe" "$packet" "$sites" > "$dest/inputs.sha256"
mpirun -np 4 "$exe" "$packet" "$sites" "$dest/transfer_mpi4" 12 10 0 0 > "$dest/transfer_mpi4/contract.json" 2> "$dest/transfer_mpi4/stderr.log"
cmp "$dest/transfer_mpi4/transfer_moments.csv" "$root/rk4_diagnostics/curved_transfer_v3/coarse_mpi4/transfer_moments.csv"
cmp "$dest/transfer_mpi4/angular_uncertainty.csv" "$root/rk4_diagnostics/curved_transfer_v3/coarse_mpi4/angular_uncertainty.csv"
printf 'AUDIT_COMPLETE\n' > "$dest/status.txt"
