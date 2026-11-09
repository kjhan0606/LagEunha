#!/usr/bin/env bash
set -euo pipefail
project=$(cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$project"
threads=${BUILD_JOBS:-2}
if [[ -n "${SLURM_JOB_ID:-}" ]];then
  cpu_limit=${SLURM_CPUS_PER_TASK:-1}
elif [[ "${BUILD_LOCAL:-0}" == 1 ]];then
  # Explicit CPU-only compilation on an authorized development host.
  cpu_limit=$(getconf _NPROCESSORS_ONLN)
else
  echo 'Use a Slurm allocation or explicitly set BUILD_LOCAL=1.' >&2;exit 2
fi
test "$threads" -ge 1 && test "$threads" -le "$cpu_limit"
cuda_arch=${CUDA_ARCH:-90}
[[ "$cuda_arch" =~ ^[0-9]+$ ]] || { echo 'Invalid CUDA_ARCH.' >&2;exit 2; }
if [[ "${BUILD_RESUME_LINK_ONLY:-0}" == 1 ]];then
  sha256sum -c ../resume_sources_and_objects.sha256 > ../checked_resume_objects.log
else
  mkdir .native-u64
fi
build="$project/.native-u64"
source_dir="$project/Exam/NativeU64"
profile=${NATIVE_PROFILE:-legacy}
material_main_defines=(); material_boundary_defines=(); material_geometry_defines=(); material_host_defines=(); material_host_includes=(); material_link_flags=(); material_objects=(); material_weight_defines=(); material_project_defines=()
case "$profile" in
 legacy) ;;
 material) source "$source_dir/material_profile.sh" ;;
 *) echo "Unknown NATIVE_PROFILE: $profile" >&2; exit 2 ;;
esac
# A link-only resume must not mix legacy geometry with material host objects.
# Historical manifests without a profile marker belong to the legacy build.
if [[ "${BUILD_RESUME_LINK_ONLY:-0}" == 1 ]];then
  previous_profile=legacy
  if [[ -f "$build/PROFILE" ]];then previous_profile=$(cat "$build/PROFILE");fi
  [[ "$profile" == "$previous_profile" ]] || {
    echo "Cannot resume $previous_profile objects as $profile; use a fresh checkout." >&2; exit 2;
  }
else
  printf '%s\n' "$profile" > "$build/PROFILE"
fi
flags=(-g -O2 -fno-strict-aliasing -qopenmp -DINTEL -DXYZDBL -DGOTPM -DUSE_CUDA -DGFS_WEIGHT_DYNAMICS -DGFS_FIXEDPOLY_SIDES=1024)
inc=(-I"$source_dir" -I"$project/Exam" -I"$project" -I"$project/Voro" -I/home/kjhan/local/include -I"$project/RkS" -I"$project/Cosmos" -I"$project/Timer" -I"$project/MpiAux" -I"$project/Colorize" -I"$project/OST" -I"$project/Exam/KH")
# No inherited objects or archives: every archive below is compiled in this checkout.
common=("CC=mpiicx" "FC=mpiifx -nofor-main" "F90C=ifx -nofor-main" "OPT=-g -O2 -fno-strict-aliasing -qopenmp" "specialrules=-DINTEL -DXYZDBL -DGOTPM -DUSE_CUDA -DGFS_WEIGHT_DYNAMICS -DGFS_FIXEDPOLY_SIDES=1024" "FFTW=/home/kjhan/local" "USE_CUDA=1" "CUDA=/opt/ohpc/pub/cuda/13.0.2")
archives=(MpiAux/libmpiaux.a Cosmos/CosmosIC/CAMB_Sources/libcamb.a Params/libheader.a Cosmos/libcosmology.a Cosmos/CosmosIC/libtwolpt.a SPH/libsph.a Voro/libmyvoro.a RkS/librks.a OST/libost.a RECFAST/librecfast.a SPH/JShin/cooling/libget_lambda.a SPH/JShin/ssp/libnew_ssp.a SPH/JShin/Extern/libextforce.a Cosmos/CosmosRW/libkjhrw.a Timer/libtimerutil.a FFT/libfft.a DenRW/libdenrw.a ParallelIO/libparallelIO.a Colorize/libcolor.a Exam/RT_LF/librtold.a Exam/RT/librt.a Exam/Kepler/libkp.a Exam/MkGlass2D/libglass2d.a Exam/Cylinder/libcyl.a Exam/Sedov2D/libsedov2d.a)
if [[ "${BUILD_RESUME_LINK_ONLY:-0}" != 1 ]];then
for lib in "${archives[@]}";do
    test ! -e "$lib" || { echo "Use a fresh checkout; existing archive $lib" >&2;exit 2; }
done
mpiicx "${flags[@]}" "${inc[@]}" "$source_dir/layout_check.c" -o "$build/layout_native.exe"
"$build/layout_native.exe" > "$build/layout_native.log"
mpiicx "${flags[@]}" -UGFS_WEIGHT_DYNAMICS "${inc[@]}" "$source_dir/layout_check.c" -o "$build/layout_legacy.exe"
"$build/layout_legacy.exe" > "$build/layout_legacy.log"
for lib in "${archives[@]}";do
    jobs=$threads
    [[ "$lib" == Cosmos/CosmosIC/CAMB_Sources/* ]] && jobs=1
    # These legacy Makefiles have a %.c generator with no prerequisites.
    # Built-in makefile-remake rules plus -B try to compile Makefile.c.
    # A fresh tree already forces object compilation; use explicit rules only.
    boundary_extra=()
    if [[ "$profile" == material && "$lib" == Exam/Cylinder/libcyl.a ]];then
      boundary_extra=("specialrules=-DINTEL -DXYZDBL -DGOTPM -DUSE_CUDA -DGFS_WEIGHT_DYNAMICS -DGFS_FIXEDPOLY_SIDES=1024 ${material_boundary_defines[*]}")
    fi
    make -r -C "${lib%/*}" -j "$jobs" "${common[@]}" "${boundary_extra[@]}" "${lib##*/}"
    test -f "$lib"
done
/opt/ohpc/pub/cuda/13.0.2/bin/nvcc -O3 --use_fast_math --threads 1 -DXYZDBL -DGOTPM -DUSE_CUDA -DGFS_WEIGHT_DYNAMICS -DGFS_VARPRESS_W0 -DGFS_FIXEDPOLY_SIDES=1024 -gencode "arch=compute_$cuda_arch,code=sm_$cuda_arch" "${material_geometry_defines[@]}" "${inc[@]}" -c "$source_dir/exam_gpu.cu" -o "$build/exam_gpu.o"
/opt/ohpc/pub/cuda/13.0.2/bin/nvcc -O3 --use_fast_math --threads 1 -DXYZDBL -DGOTPM -DUSE_CUDA -DGFS_WEIGHT_DYNAMICS -DGFS_VARPRESS_W0 -DGFS_FIXEDPOLY_SIDES=1024 -gencode "arch=compute_$cuda_arch,code=sm_$cuda_arch" "${inc[@]}" -c "$source_dir/gfs_unified_gpu.cu" -o "$build/gfs_unified_gpu.o"
mpiicx "${flags[@]}" "${inc[@]}" -c "$project/Exam/gfs_viscosity_implicit.c" -o "$build/gfs_viscosity_implicit.o"
for unit in exam exam_gpu_extract mpirks.exam2d exam3d_gfs gfs_circle_observer gfs_circle_dense kh util gfs_weight_checkpoint;do
    extra=(); extra_inc=()
    if [[ "$profile" == material && "$unit" == exam ]];then extra=("${material_main_defines[@]}");fi
    if [[ "$profile" == material && "$unit" == exam_gpu_extract ]];then extra=("${material_host_defines[@]}"); extra_inc=("${material_host_includes[@]}");fi
    mpiicx "${flags[@]}" "${extra[@]}" "${extra_inc[@]}" "${inc[@]}" -c "$source_dir/$unit.c" -o "$build/$unit.o"
done
fi
if [[ "$profile" == material ]];then
 material_objects=("$build/material_cuda_rhs.o" "$build/material_cuda_weight.o" "$build/material_cuda_project.o")
 if [[ "${BUILD_RESUME_LINK_ONLY:-0}" != 1 ]];then
  /opt/ohpc/pub/cuda/13.0.2/bin/nvcc -O2 -arch=sm_80 "${material_rhs_defines[@]}" -I"$project/Exam" -c "$source_dir/gfs_material_cuda_rhs.cu" -o "$build/material_cuda_rhs.o"
  /opt/ohpc/pub/cuda/13.0.2/bin/nvcc -O2 -arch=sm_80 "${material_weight_defines[@]}" -I"$project/Exam" -c "$source_dir/gfs_material_cuda_weight.cu" -o "$build/material_cuda_weight.o"
  /opt/ohpc/pub/cuda/13.0.2/bin/nvcc -O2 --fmad=false -arch=sm_80 "${material_project_defines[@]}" -I"$project/Exam" -c "$source_dir/gfs_material_cuda_project.cu" -o "$build/material_cuda_project.o"
 fi
fi
mpiicx "${flags[@]}" "${inc[@]}" -I"$project/FFT" -I"$project/Params" -I"$project/SPH" -I"$project/ParallelIO" -I"$project/Exam/RT" -I"$project/Exam/Kepler" -I"$project/Exam/MkGlass2D" -I"$project/Exam/Cylinder" -I"$project/Exam/Sedov2D" -c eunha2.c -o "$build/eunha2.o"
ar rcs "$build/libexam.a" "$build/exam.o" "$build/exam_gpu.o" "$build/gfs_unified_gpu.o" "$build/gfs_viscosity_implicit.o" "$build/exam_gpu_extract.o" "$build/mpirks.exam2d.o" "$build/exam3d_gfs.o" "$build/gfs_circle_observer.o" "$build/gfs_circle_dense.o" "${material_objects[@]}"
ar rcs "$build/libkh.a" "$build/kh.o" "$build/util.o" "$build/gfs_weight_checkpoint.o"
# Preserve the qualified historical archive order; use only fresh local libraries.
libpaths=(-L"$build")
for lib in "${archives[@]}";do libpaths+=(-L"$project/${lib%/*}");done
libs=(-lheader -lcosmology -ltwolpt -lkjhrw -lrks -ltimerutil -lsph -lrecfast -ldenrw -lmpiaux -lparallelIO -lkh -lrt -lrtold -lkp -lglass2d -lcyl -lsedov2d -lexam -lost -lcolor -lmyvoro -lget_lambda -lnew_ssp -lextforce -lfft -lcamb)
mpiifx -nofor-main -g -O2 -qopenmp -DINTEL -DXYZDBL -DGOTPM "$build/eunha2.o" "${libpaths[@]}" "${libs[@]}" -L/home/kjhan/local/lib -lfftw3f_mpi -lfftw3f -lfftw3f_threads -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lcudadevrt -lstdc++ "${material_link_flags[@]}" -o "$build/eunha2_native_u64.exe"
sha256sum "$build/eunha2_native_u64.exe" "$build/"*.o "$build/"*.a "${archives[@]}" > "$build/SHA256SUMS"
printf 'Full native U64 source build complete. Physical parity not yet qualified.\n' > "$build/COMPILE_COMPLETE"
