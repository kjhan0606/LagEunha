#!/usr/bin/env bash
#SBATCH --job-name=power-edge-v2-sm90
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/material_compact_actual_source_v1/slurm-power-edge-%j.out
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
test "$#" = 2
first=$1
last=$2
[[ "$first" =~ ^[0-9]+$ && "$last" =~ ^[0-9]+$ ]]
(( last > first && last - first <= 50000 ))
run=$dest/power_edge_chunk_${first}_${last}_h200sm90_v2b
test -n "${SLURM_JOB_ID:-}"
source "$transition/run_env.sh"
export I_MPI_CXX=g++ I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0
export OMP_NUM_THREADS=4 OMP_PROC_BIND=false
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
cd "$dest"
test "$(cat context_native_exit_code.txt)" = 0
mkdir "$run"
sha256sum -c context_geometry.sha256 > "$run/cached_geometry_check.txt"
sha256sum -c context_source_input.sha256 > "$run/context_input_check.txt"
sha256sum -c lageunha_actual_source_input.sha256 > "$run/accepted_input_check.txt"
sha256sum -c "$dest/weight_direction_v2/weights.sha256" > "$run/weight_check.txt"
sha256sum probe_power_parallel_v2.cpp run_power_v49_edge_chunk_h200_sm90_v2b.sbatch \
  "$root/material_compact_rezone_parallel_v2.cpp" \
  "$root/material_compact_chart_topology.h" \
  "$root/material_compact_chart_jacobian.h" \
  "$root/material_compact_chart_area.h" \
  "$root/material_compact_power_cell.h" \
  "$root/material_compact_power_curve.h" \
  "$root/material_compact_overlap_boundary_v2.h" \
  "$transition/Exam/libexam.a" \
  "$transition/Exam/material_full_initial_cuda.cu" \
  "$dest/weight_direction_v2/weights.csv" > "$run/source.sha256"
mkdir "$run/Exam"
cp "$transition/Exam/libexam.a" "$run/Exam/libexam.a"
includes=(-I"$transition/Exam" -I"$original" -I"$original/Exam" -I"$original/Voro" \
  -I/home/kjhan/local/include -I"$original/RkS" -I"$original/Cosmos" \
  -I"$original/Timer" -I"$original/MpiAux" -I"$original/Colorize" -I"$original/OST")
nvcc -O2 -std=c++17 --fmad=false -DXYZDBL -DGOTPM -DUSE_CUDA \
  -gencode arch=compute_90,code=sm_90 "${includes[@]}" \
  -c "$transition/Exam/material_full_initial_cuda.cu" \
  -o "$run/Exam/material_full_initial_cuda.o" > "$run/build_sm90.log" 2>&1
ar rcs "$run/Exam/libexam.a" "$run/Exam/material_full_initial_cuda.o"
sha256sum "$run/Exam/material_full_initial_cuda.o" "$run/Exam/libexam.a" > "$run/h200_cuda_binary.sha256"
mpicxx -std=gnu++17 -fext-numeric-literals -O2 -ffunction-sections -fdata-sections \
  -fopenmp -fno-fast-math -ffp-contract=off -DXYZDBL -DGOTPM -DUSE_CUDA \
  -I "$original/Exam" -I /opt/ohpc/pub/cuda/13.0.2/include \
  probe_power_parallel_v2.cpp -Wl,--gc-sections -L "$run/Exam" -lexam \
  -L /opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lquadmath \
  -o "$run/probe_power"
sha256sum "$run/probe_power" > "$run/binary.sha256"
TIMEFORMAT='elapsed_seconds=%3R'
set +e
time timeout --signal=TERM --kill-after=60s 43200s \
  mpirun -np 4 "$run/probe_power" "$initial" \
  "$transition/immutable_payload.sha256" "$producer" "$whole" \
  "$reference" "$accepted" "$dest/context_v1" \
  "$dest/weight_direction_v2/weights.csv" 1 boundary_parallel_range \
  "$first" "$last" \
  </dev/null > "$run/native.log" 2>&1
result=$?
set -e
printf '%s\n' "$result" > "$run/native_exit_code.txt"
sha256sum -c "$run/source.sha256" > "$run/final_source_check.txt"
if test "$result" = 0; then
  sha256sum "$run"/source_edge_wires_rank{0,1,2,3}.csv \
    "$run"/source_edge_manifest_rank{0,1,2,3}.csv > "$run/source_edge_wires.sha256"
fi
tail -35 "$run/native.log"
exit "$result"
