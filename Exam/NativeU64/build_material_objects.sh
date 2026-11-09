#!/usr/bin/env bash
# Compile the material modules from this checkout, into a fresh output directory.
set -euo pipefail
project=$(cd -- "$(dirname -- "$0")/../.." && pwd)
source_dir="$project/Exam/NativeU64"
build=${1:?Usage: build_material_objects.sh absolute-output-directory}
[[ "$build" = /* ]] || { echo 'Use an absolute output directory.' >&2; exit 2; }
mkdir -- "$build"
source "$source_dir/material_profile.sh"
cuda_home=${CUDA_HOME:-/opt/ohpc/pub/cuda/13.0.2}
cuda_arch=${CUDA_ARCH:-90}
[[ "$cuda_arch" =~ ^[0-9]+$ ]] || exit 2
flags=(-g -O2 -fno-strict-aliasing -qopenmp -DINTEL -DXYZDBL -DGOTPM -DUSE_CUDA -DGFS_WEIGHT_DYNAMICS -DGFS_FIXEDPOLY_SIDES=1024)
inc=(-I"$source_dir" -I"$project/Exam" -I"$project" -I"$project/Voro" -I/home/kjhan/local/include -I"$project/RkS" -I"$project/Cosmos" -I"$project/Timer" -I"$project/MpiAux" -I"$project/Colorize" -I"$project/OST" -I"$project/Exam/KH")
{
 set -x
 mpiicx "${flags[@]}" "${inc[@]}" -c "$source_dir/exam.c" -o "$build/exam.o"
 mpiicx "${flags[@]}" "${material_host_defines[@]}" "${material_host_includes[@]}" "${inc[@]}" -c "$source_dir/exam_gpu_extract.c" -o "$build/exam_gpu_extract.o"
 "$cuda_home/bin/nvcc" -O3 --use_fast_math --threads 1 -gencode "arch=compute_$cuda_arch,code=sm_$cuda_arch" -DXYZDBL -DGOTPM -DUSE_CUDA -DGFS_WEIGHT_DYNAMICS -DGFS_VARPRESS_W0 -DGFS_FIXEDPOLY_SIDES=1024 "${material_geometry_defines[@]}" "${inc[@]}" -c "$source_dir/exam_gpu.cu" -o "$build/exam_gpu.o"
 # Preserve the validated face-response arithmetic and device target during
 # this source-port gate. Architecture retuning needs a separate comparison.
 "$cuda_home/bin/nvcc" -O2 -arch=sm_80 "${material_rhs_defines[@]}" -I"$project/Exam" -c "$source_dir/gfs_material_cuda_rhs.cu" -o "$build/material_cuda_rhs.o"
 "$cuda_home/bin/nvcc" -O2 -arch=sm_80 "${material_weight_defines[@]}" -I"$project/Exam" -c "$source_dir/gfs_material_cuda_weight.cu" -o "$build/material_cuda_weight.o"
 "$cuda_home/bin/nvcc" -O2 --fmad=false -arch=sm_80 "${material_project_defines[@]}" -I"$project/Exam" -c "$source_dir/gfs_material_cuda_project.cu" -o "$build/material_cuda_project.o"
 ar rcs "$build/libmaterial_native.a" "$build/exam.o" "$build/exam_gpu_extract.o" "$build/exam_gpu.o" "$build/material_cuda_rhs.o" "$build/material_cuda_weight.o" "$build/material_cuda_project.o"
 set +x
} > "$build/build.log" 2>&1
sha256sum "$build/"*.o "$build/"*.a > "$build/SHA256SUMS"
printf 'Material native objects compiled from the current source tree.\n'
