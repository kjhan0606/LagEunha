#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
audit=$root/rk4_diagnostics/dg_projection_${version}_build
mkdir "$audit"
source "$root/run_env.sh"
export OMP_NUM_THREADS=1
cd /gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src
sha256sum Exam/exam.c Exam/exam_gpu.cu Exam/exam_gpu.h Exam/exam_gpu_extract.c Exam/KH/kh.c Exam/KH/Makefile Exam/gfs_ssprk3.inc Exam/gfs_ale_flux.h Exam/gfs_ale_shadow.inc Exam/gfs_moment.h Exam/gfs_moment_density.h Exam/gfs_moment_initial.h Exam/gfs_moment_gpu.inc Exam/gfs_moment_native.inc Exam/gfs_dg_euler.h Exam/gfs_dg_gpu.inc Exam/gfs_dg_native.inc Exam/Makefile > "$audit/source.sha256"
cp Exam/gfs_moment.h Exam/gfs_moment_density.h Exam/gfs_moment_initial.h Exam/gfs_dg_euler.h Exam/gfs_dg_gpu.inc Exam/gfs_dg_native.inc "$audit/"
cp "$root/rk4_diagnostics/test_dg_euler_projection_v1.cpp" "$audit/test_dg_euler_projection.cpp"
g++ -O2 -std=c++14 -I"$audit" "$audit/test_dg_euler_projection.cpp" -o "$audit/test_dg_euler_projection.exe"
"$audit/test_dg_euler_projection.exe" > "$audit/contract.json"
make -C Exam -j1 USE_CUDA=1 libexam.a > "$audit/build_exam.log" 2>&1
make -C Exam/KH -j1 USE_CUDA=1 libkh.a > "$audit/build_kh.log" 2>&1
make -B -j1 USE_CUDA=1 eunha2.exe > "$audit/build_link.log" 2>&1
cp eunha2.exe "$audit/eunha2_dg_projection.exe"
sha256sum "$audit/eunha2_dg_projection.exe" > "$audit/binary.sha256"
sha256sum -c "$audit/source.sha256" > "$audit/source_unchanged.txt"
printf 'BUILD_COMPLETE\n' > "$audit/build_status.txt"
