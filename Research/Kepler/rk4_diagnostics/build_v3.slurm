#!/usr/bin/env bash
#SBATCH --job-name=KepRK4FilterBuild
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/build-v3-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/build-v3-%j.err
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
audit=$root/rk4_diagnostics/v3
mkdir "$audit"
source "$root/run_env.sh"
cd /gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src
cp Exam/exam.c "$audit/exam.after.c"
cp Exam/KH/kh.c "$audit/kh.after.c"
cp Exam/exam2d.h "$audit/exam2d.after.h"
sha256sum Exam/exam.c Exam/KH/kh.c Exam/exam2d.h Exam/exam_gpu.cu Exam/gfs_velocity_limiter.h > "$audit/source.sha256"
make -C Exam -j1 USE_CUDA=1 libexam.a > "$audit/build_exam.log" 2>&1
make -C Exam/KH -j1 USE_CUDA=1 libkh.a > "$audit/build_kh.log" 2>&1
make -B -j1 USE_CUDA=1 eunha2.exe > "$audit/build_link.log" 2>&1
cp eunha2.exe "$audit/eunha2_rk4_diagnostics.exe"
sha256sum "$audit/eunha2_rk4_diagnostics.exe" > "$audit/binary.sha256"
printf 'BUILD_COMPLETE\n' > "$audit/build_status.txt"
