#!/usr/bin/env bash
#SBATCH --job-name=KepCoreRefine
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=12G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/core-native-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/core-native-%j.err
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
source "$root/run_env.sh"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 "$root/rk4_diagnostics/run_preflight_and_capture.py" --run-dir validation_v3 --binary v3/eunha2_rk4_diagnostics.exe --preflight-only
python3 "$root/run_core_refinement.py"
