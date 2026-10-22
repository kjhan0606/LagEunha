#!/usr/bin/env bash
#SBATCH --job-name=KepCoreRefineV2
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=12G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/core-native-v2-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/core-native-v2-%j.err
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
source "$root/run_env.sh"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 "$root/run_core_refinement.py" --run-dir native_v2 --binary raw_vertex_dedup/eunha2_raw_vertex_dedup.exe --preflight-only
