#!/usr/bin/env bash
#SBATCH --job-name=KepCoreIC
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/core-ic-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/core-ic-%j.err
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
source "$root/run_env.sh"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 "$root/prepare_core_refinement.py"
