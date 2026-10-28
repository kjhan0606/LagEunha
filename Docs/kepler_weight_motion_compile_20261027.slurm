#!/usr/bin/env bash
#SBATCH --job-name=LagQ2WeightCompile
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:0
#SBATCH --mem=4G
#SBATCH --time=00:05:00
set -euo pipefail
module purge
module load gnu13/13.2.0 cuda/13.0.2
sha256sum -c inputs.sha256
nvcc -std=c++14 -O2 -arch=sm_80 test.cu -o "${SLURM_TMPDIR:-/tmp}/weight_motion_test_${SLURM_JOB_ID}"
echo 'PASS: CUDA translation and link; GPU execution remains a separate test'
