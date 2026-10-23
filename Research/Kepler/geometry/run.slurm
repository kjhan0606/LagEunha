#!/bin/bash
#SBATCH --job-name=compact-actual-step5-source
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/material_compact_actual_source_v1/slurm-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/material_compact_actual_source_v1/slurm-%j.err

set -euo pipefail
module purge
module load gnu13/13.2.0 intel/mpi/2021.17 cuda/13.0.2
export I_MPI_CXX=g++ I_MPI_HYDRA_BOOTSTRAP=slurm OMP_NUM_THREADS=2 OMP_PROC_BIND=false
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
dest=$root/material_compact_actual_source_v1
original=/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src
initial=$root/material_full_initial_v2_h200_mpi4/run
producer=$root/material_full_positive_rule_v1_ada_mpi4
transition=$root/material_full_positive_transition_v2_ada_mpi4
whole=$root/material_positive_rule_whole_audit_v1_mpi4
reference=$root/material_full_reference_v2_h200_mpi4/run/material_full_reference.json
accepted=$transition/run/continuous_step5/accepted_r

cd "$dest"
sha256sum -c "$transition/immutable_payload.sha256" > immutable_input_check.txt
sha256sum -c "$transition/positive_rule_input.sha256" > positive_input_check.txt
sha256sum "$dest/probe.cpp" "$root/material_compact_rezone.cpp" \
  "$reference" "$whole/contract.json" \
  "$transition/run/continuous_step5/step.json" \
  "${accepted}0.bin" "${accepted}1.bin" "${accepted}2.bin" "${accepted}3.bin" \
  > actual_source_input.sha256
mpicxx -std=gnu++17 -fext-numeric-literals -O2 -ffunction-sections -fdata-sections \
  -fopenmp -fno-fast-math -ffp-contract=off -DXYZDBL -DGOTPM -DUSE_CUDA \
  -I "$original/Exam" -I /opt/ohpc/pub/cuda/13.0.2/include \
  probe.cpp -Wl,--gc-sections -L /opt/ohpc/pub/cuda/13.0.2/lib64 \
  -lquadmath -lcudart -o probe_optimized
sha256sum probe_optimized > probe_binary.sha256
TIMEFORMAT='elapsed_seconds=%3R'
time mpirun -n 4 ./probe_optimized "$initial" \
  "$transition/immutable_payload.sha256" "$producer" "$whole" \
  "$reference" "$accepted" </dev/null
sha256sum -c actual_source_input.sha256 > final_source_input_check.txt
