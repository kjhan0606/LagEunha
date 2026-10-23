#!/usr/bin/env bash
#SBATCH --job-name=power-triangle-route-root60
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/material_compact_actual_source_v1/slurm-power-triangle-route-%j.out
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
dest=$root/material_compact_actual_source_v1
original=/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src
initial=$root/material_full_initial_v2_h200_mpi4/run
producer=$root/material_full_positive_rule_v1_ada_mpi4
transition=$root/material_full_positive_transition_v2_ada_mpi4
whole=$root/material_positive_rule_whole_audit_v1_mpi4
reference=$root/material_full_reference_v2_h200_mpi4/run/material_full_reference.json
accepted=$transition/run/continuous_step5/accepted_r
test "$#" = 3
first=$1
last=$2
chunk=$3
[[ "$first" =~ ^[0-9]+$ && "$last" =~ ^[0-9]+$ ]]
(( last > first && last - first <= 50000 ))
test -n "${SLURM_JOB_ID:-}"
cd "$dest"
test "$(cat "$chunk/native_exit_code.txt")" = 0
canonical=$dest/power_edge_reconcile_${first}_${last}_root60_v1
test "$(cat "$canonical/native_exit_code.txt")" = 0
run=$dest/power_edge_triangle_route_chunk_${first}_${last}_root60_v1
h200lib=$dest/power_edge_chunk_1000000_1050000_h200sm90_root60_omp8/Exam/libexam.a
source "$transition/run_env.sh"
export I_MPI_CXX=g++ I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0
export OMP_NUM_THREADS=4 OMP_PROC_BIND=false
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
cd "$dest"
sha256sum -c "$chunk/source_edge_wires.sha256" >/dev/null
sha256sum -c "$canonical/canonical.sha256" >/dev/null
sha256sum -c "$dest/weight_direction_v2/weights.sha256" >/dev/null
sha256sum -c "$dest/power_edge_chunk_1000000_1050000_h200sm90_root60_omp8/h200_cuda_binary.sha256" >/dev/null
mkdir "$run"
sha256sum probe_power_triangle_route_root60_v1.cpp run_power_triangle_route_root60_h200_v1.sbatch \
  "$root/material_compact_rezone_parallel_v2.cpp" \
  "$root/material_compact_overlap_boundary_v2.h" \
  "$root/material_compact_overlap_cuts.h" \
  "$root/material_compact_overlap_triangle_route.h" \
  "$root/material_compact_overlap_events.h" \
  "$h200lib" \
  "$dest/weight_direction_v2/weights.csv" \
  "$chunk"/source_edge_manifest_rank{0,1,2,3}.csv \
  "$canonical"/canonical_owner{0,1,2,3}.csv > "$run/source.sha256"
mpicxx -std=gnu++17 -fext-numeric-literals -O2 -ffunction-sections -fdata-sections \
  -fopenmp -fno-fast-math -ffp-contract=off -DXYZDBL -DGOTPM -DUSE_CUDA \
  -I "$original/Exam" -I /opt/ohpc/pub/cuda/13.0.2/include \
  probe_power_triangle_route_root60_v1.cpp -Wl,--gc-sections "$h200lib" \
  -L /opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lquadmath \
  -o "$run/probe_power_triangle_route"
sha256sum "$run/probe_power_triangle_route" > "$run/binary.sha256"
TIMEFORMAT='elapsed_seconds=%3R'
set +e
time timeout --signal=TERM --kill-after=60s 7200s \
  mpirun -np 4 "$run/probe_power_triangle_route" "$initial" \
  "$transition/immutable_payload.sha256" "$producer" "$whole" \
  "$reference" "$accepted" "$dest/context_v1" \
  "$dest/weight_direction_v2/weights.csv" 1 boundary_triangle_route_range \
  "$first" "$last" "$chunk" \
  </dev/null > "$run/native.log" 2>&1
result=$?
set -e
printf '%s\n' "$result" > "$run/native_exit_code.txt"
sha256sum -c "$run/source.sha256" > "$run/final_source_check.txt"
tail -25 "$run/native.log"
exit "$result"
