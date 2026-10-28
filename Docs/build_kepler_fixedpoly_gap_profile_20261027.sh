#!/usr/bin/env bash
set -euo pipefail

source_dir=/tmp/Kepler_sparse_fixedpoly_profile_20261027
baseline_dir=/tmp/Kepler_sparse_fixedpoly_20261027
expected_gpu_source=2194c4aa81ea0eee030c5b07955794d0b7dc8d25dbc32f5049d14c71add67865
expected_gpu_object=e0d30024a9afc5abab3965ba39d01d04ea5171b01d527b7f9b905c5afd2ccd6c
expected_profile_patch=Docs/KEPLER_RK4_GAP_PROFILE_2026-10-27.patch

test "$(hostname)" = syntax
test -f "$source_dir/Exam/exam.c"
test -f "$source_dir/Exam/exam_gpu.cu"
test -f "$expected_profile_patch"
actual_gpu_source=$(sha256sum "$source_dir/Exam/exam_gpu.cu")
test "${actual_gpu_source%% *}" = "$expected_gpu_source"
actual_gpu_object=$(sha256sum "$baseline_dir/Exam/exam_gpu.o")
test "${actual_gpu_object%% *}" = "$expected_gpu_object"

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2

make -C "$source_dir/Exam" -j2 USE_CUDA=1 GFS_VARPRESS_W0=1 libexam.a
# Only exam.c carries the timer changes. Reuse the *identical* accepted
# GPU object even if make rebuilt it while refreshing libexam.a.
cp "$baseline_dir/Exam/exam_gpu.o" "$source_dir/Exam/exam_gpu.o"
ar rcs "$source_dir/Exam/libexam.a" "$source_dir/Exam/exam_gpu.o"
make -C "$source_dir" -j2 USE_CUDA=1 GFS_VARPRESS_W0=1 this

sha256sum "$source_dir/Exam/exam.c" "$source_dir/Exam/exam_gpu.o" "$source_dir/eunha2.exe"
