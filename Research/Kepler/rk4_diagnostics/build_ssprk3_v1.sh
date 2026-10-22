#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
audit=$root/rk4_diagnostics/ssprk3_v1_build
mkdir "$audit"
source "$root/run_env.sh"
export OMP_NUM_THREADS=1
cd /gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src
sha256sum Exam/exam.c Exam/gfs_ssprk3.inc Exam/exam_gpu.cu Exam/exam_gpu.h Exam/exam_gpu_extract.c Exam/KH/kh.c > "$audit/source.sha256"
cp Exam/gfs_ssprk3.inc "$audit/gfs_ssprk3.inc"
make -C Exam -j1 USE_CUDA=1 libexam.a > "$audit/build_exam.log" 2>&1
make -C Exam/KH -j1 USE_CUDA=1 libkh.a > "$audit/build_kh.log" 2>&1
make -B -j1 USE_CUDA=1 eunha2.exe > "$audit/build_link.log" 2>&1
cp eunha2.exe "$audit/eunha2_ssprk3.exe"
sha256sum "$audit/eunha2_ssprk3.exe" > "$audit/binary.sha256"
sha256sum -c "$audit/source.sha256" > "$audit/source_unchanged.txt"
printf 'BUILD_COMPLETE\n' > "$audit/build_status.txt"
