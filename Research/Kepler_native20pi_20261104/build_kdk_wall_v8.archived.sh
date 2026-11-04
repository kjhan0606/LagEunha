#!/usr/bin/env bash
set -euo pipefail
pkg=$(cd -- "$(dirname -- "$0")" && pwd)
root=/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029
base="$root/native_u64_w0_long_20261101_v2"
state="$root/native_weight_state_20261101_v1"
prior="$root/circle_boundary_event_integration_20261101_v1"
parent="$root/polygon_refinement_native_t35_20261101_v1"
project="$state/project"

test -f "$state/STATE_IO_COMPLETE"
(cd "$base"; sha256sum -c SHA256SUMS) > "$pkg/checked_base_source.log"
(cd "$base"; sha256sum -c build/SHA256SUMS) > "$pkg/checked_base_build.log"
(cd "$state"; sha256sum -c SHA256SUMS) > "$pkg/checked_state_source.log"
(cd "$state"; sha256sum -c build/SHA256SUMS) > "$pkg/checked_state_build.log"
(cd "$pkg"; sha256sum -c source_manifest.sha256) > "$pkg/checked_candidate_source.log"

build="$pkg/build_kdk_wall_v8"
mkdir "$build"
for item in exam_gpu.o exam_gpu_extract.o libexam.a libkh.a; do
    cp "$base/build/$item" "$build/$item"
done
cuda_ref="/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/native_u64_kdk_individual_20261104_v1/kdk_block2_pilot_v1/build_rerun_manual_a10b"
current_cuda_source=$(sha256sum "$pkg/source_overlay/exam_gpu.cu" | awk '{print $1}')
if [[ -f "$cuda_ref/exam_gpu.source.sha256" &&
      "$(cat "$cuda_ref/exam_gpu.source.sha256")" == "$current_cuda_source" ]]; then
    cp "$cuda_ref/exam_gpu.o" "$build/exam_gpu.o"
    printf '%s\n' "$current_cuda_source" > "$build/exam_gpu.source.sha256"
else
/opt/ohpc/pub/cuda/13.0.2/bin/nvcc -O3 --use_fast_math --threads 1 \
  -DXYZDBL -DGOTPM -DUSE_CUDA -DGFS_WEIGHT_DYNAMICS -DGFS_VARPRESS_W0 \
  -DGFS_FIXEDPOLY_SIDES=1024 \
  -gencode arch=compute_86,code=sm_86 \
  -gencode arch=compute_89,code=sm_89 \
  -gencode arch=compute_90,code=sm_90 \
  -I"$pkg/source_overlay" -I"$prior/source_overlay" -I"$parent/source_overlay" \
  -I"$project/Exam" -I"$project" -I"$project/Voro" -I/home/kjhan/local/include \
  -I"$project/RkS" -I"$project/Cosmos" -I"$project/Timer" -I"$project/MpiAux" \
  -I"$project/Colorize" -I"$project/OST" -I"$project/Exam/KH" \
  -c "$pkg/source_overlay/exam_gpu.cu" -o "$build/exam_gpu.o"
  printf '%s\n' "$current_cuda_source" > "$build/exam_gpu.source.sha256"
fi
mpiicx -g -O2 -fno-strict-aliasing -qopenmp -DINTEL -DXYZDBL -DGOTPM -DUSE_CUDA -DGFS_WEIGHT_DYNAMICS -DGFS_FIXEDPOLY_SIDES=1024 \
  -I"$pkg/source_overlay" -I"$prior/source_overlay" -I"$parent/source_overlay" \
  -I"$project/Exam" -I"$project" -I"$project/Voro" -I/home/kjhan/local/include \
  -I"$project/RkS" -I"$project/Cosmos" -I"$project/Timer" -I"$project/MpiAux" \
  -I"$project/Colorize" -I"$project/OST" -I"$project/Exam/KH" \
  -c "$pkg/source_overlay/exam.c" -o "$build/exam.o"
mpiicx -g -O2 -fno-strict-aliasing -qopenmp -DINTEL -DXYZDBL -DGOTPM -DUSE_CUDA -DGFS_WEIGHT_DYNAMICS -DGFS_FIXEDPOLY_SIDES=1024 \
  -I"$pkg/source_overlay" -I"$prior/source_overlay" -I"$parent/source_overlay" \
  -I"$project/Exam" -I"$project" -I"$project/Voro" -I/home/kjhan/local/include \
  -I"$project/RkS" -I"$project/Cosmos" -I"$project/Timer" -I"$project/MpiAux" \
  -I"$project/Colorize" -I"$project/OST" -I"$project/Exam/KH" \
  -c "$pkg/source_overlay/kh.c" -o "$build/kh.o"
(cd "$build"; ar r libexam.a exam.o exam_gpu.o; ranlib libexam.a; ar r libkh.a kh.o; ranlib libkh.a)
libs="-lheader -lcosmology -ltwolpt -lkjhrw -lrks -ltimerutil -lsph -lrecfast -ldenrw -lmpiaux -lparallelIO -lkh -lrt -lrtold -lkp -lglass2d -lcyl -lsedov2d"
tail_libs="-lost -lcolor -lmyvoro -lget_lambda -lnew_ssp -lextforce -lfft -lcamb -L/home/kjhan/local/lib -lfftw3f_mpi -lfftw3f -lfftw3f_threads -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lcudadevrt -lstdc++"
mpiifx -nofor-main -g -O2 -qopenmp -DINTEL -DXYZDBL -DGOTPM \
  "$parent/link_inputs/eunha2.o" -L"$build" -L"$parent/link_inputs" $libs \
  -lexam $tail_libs -o "$build/weight_native_U64_KDK.exe"
strings "$build/weight_native_U64_KDK.exe" > "$build/binary_strings.txt"
grep -q 'U64_KDK_GUARD' "$build/binary_strings.txt"
grep -q 'U64_KDK_AUDIT' "$build/binary_strings.txt"
grep -q 'KDK_LTS_BLOCK2' "$build/binary_strings.txt"
sha256sum "$build/weight_native_U64_KDK.exe" "$build/exam.o" "$build/kh.o" \
  "$build/exam_gpu.o" "$build/exam_gpu_extract.o" > "$build/SHA256SUMS"
printf 'KDK fixed-polygon event candidate compiled; fatbin sm86/sm89/sm90; physics pilot pending.\n' > "$pkg/COMPILE_COMPLETE_v8"
