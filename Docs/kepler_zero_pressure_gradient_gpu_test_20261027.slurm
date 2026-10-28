#!/usr/bin/env bash
#SBATCH --job-name=LagQ2ZeroPGPU
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
test_dir=${SLURM_TMPDIR:-/tmp}/kepler_zero_gradient_gpu_${SLURM_JOB_ID}
mkdir -p "$test_dir/old" "$test_dir/fixed"
trap 'rm -rf "$test_dir"' EXIT
for sub in old fixed; do
  tar -xf source.tar -C "$test_dir/$sub" \
    ./Exam/gfs_energy_neutral_probe.cu ./Exam/kepler_gpu_jacobi_probe_20261027.cu \
    ./Exam/exam_gpu.h ./Exam/Makefile ./Exam/gfs_moment_initial.h \
    ./Exam/gfs_moment_density.h ./Exam/gfs_moment.h ./Exam/gfs_dg_euler.h \
    ./Exam/gfs_enriched_euler.h ./Exam/gfs_enriched_coupling.h
done
python3 patch_zero.py "$test_dir/fixed" > patch_receipt.json
for sub in old fixed; do
  nvcc -std=c++14 -O2 --use_fast_math -arch=sm_86 -DUSE_CUDA -I "$test_dir/$sub/Exam" \
    test_kepler_zero_pressure_gradient_cuda_20261027.cu -o "$test_dir/$sub/test_zero"
done
"$test_dir/old/test_zero" old | tee old_result.log
"$test_dir/fixed/test_zero" fixed | tee fixed_result.log
