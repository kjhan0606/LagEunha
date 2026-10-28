#!/usr/bin/env bash
#SBATCH --job-name=LagPolyMomentGPU
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:A10:1
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH --nice=0
set -euo pipefail
module purge
module load gnu13/13.2.0 cuda/13.0.2
sha256sum -c inputs.sha256
nvidia-smi -L
nvcc -std=c++14 -O2 -arch=sm_86 -I . test_gfs_polygon_moment_cuda_20261027.cu -o polygon_moment_test
./polygon_moment_test | tee result.log
