#!/usr/bin/env bash
#SBATCH --job-name=KepV4Regression
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=12G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/native-v4-regression-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/native-v4-regression-%j.err
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
source "$root/run_env.sh"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 "$root/core_mpi_preflight.py"
python3 "$root/rk4_diagnostics/replay_failure.py" --run-dir failure_geometry_v4 --binary raw_vertex_dedup/v4/eunha2_raw_vertex_dedup.exe
