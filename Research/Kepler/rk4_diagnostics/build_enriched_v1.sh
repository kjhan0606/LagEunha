#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
version=${1:-v1}
audit=$root/rk4_diagnostics/enriched_${version}_build
mkdir "$audit"
source "$root/run_env.sh"
export OMP_NUM_THREADS=1
cd /gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src
sha256sum Exam/exam.c Exam/exam_gpu.cu Exam/exam_gpu.h Exam/exam_gpu_extract.c Exam/KH/kh.c Exam/KH/Makefile Exam/gfs_ssprk3.inc Exam/gfs_ale_flux.h Exam/gfs_ale_fan.h Exam/gfs_ale_shadow.inc Exam/gfs_moment.h Exam/gfs_moment_density.h Exam/gfs_moment_initial.h Exam/gfs_moment_gpu.inc Exam/gfs_moment_native.inc Exam/gfs_dg_euler.h Exam/gfs_dg_gpu.inc Exam/gfs_dg_native.inc Exam/gfs_enriched_euler.h Exam/gfs_enriched_gpu.inc Exam/gfs_enriched_native.inc Exam/gfs_enriched_coupling.h Exam/gfs_enriched_coupling_native.inc Exam/gfs_enriched_coupling_gpu.inc Exam/Makefile > "$audit/source.sha256"
cp Exam/gfs_moment.h Exam/gfs_moment_density.h Exam/gfs_moment_initial.h Exam/gfs_dg_euler.h Exam/gfs_ale_flux.h Exam/gfs_ale_fan.h Exam/gfs_enriched_euler.h Exam/gfs_enriched_gpu.inc Exam/gfs_enriched_native.inc Exam/gfs_enriched_coupling.h Exam/gfs_enriched_coupling_native.inc Exam/gfs_enriched_coupling_gpu.inc "$audit/"
cp "$root/rk4_diagnostics/test_enriched_euler_v1.cpp" "$audit/"
cp "$root/rk4_diagnostics/test_weight_fan_v1.cpp" "$audit/"
g++ -O2 -std=c++14 -I"$audit" "$audit/test_weight_fan_v1.cpp" -o "$audit/test_weight_fan.exe" > "$audit/build_fan.log" 2>&1
"$audit/test_weight_fan.exe" > "$audit/fan_contract.json" 2> "$audit/fan_stderr.log"
g++ -O2 -std=c++14 -I"$audit" "$audit/test_enriched_euler_v1.cpp" -o "$audit/test_enriched_euler.exe" > "$audit/build_contract.log" 2>&1
"$audit/test_enriched_euler.exe" "$root/rk4_diagnostics/dg_projection_v2_p2/fine_fresh" > "$audit/contract.json" 2> "$audit/contract_stderr.log"
make -C Exam -j1 USE_CUDA=1 libexam.a > "$audit/build_exam.log" 2>&1
make -C Exam/KH -j1 USE_CUDA=1 libkh.a > "$audit/build_kh.log" 2>&1
make -B -j1 USE_CUDA=1 eunha2.exe > "$audit/build_link.log" 2>&1
cp eunha2.exe "$audit/eunha2_enriched.exe"
sha256sum "$audit/eunha2_enriched.exe" > "$audit/binary.sha256"
sha256sum -c "$audit/source.sha256" > "$audit/source_unchanged.txt"
printf 'BUILD_COMPLETE\n' > "$audit/build_status.txt"
