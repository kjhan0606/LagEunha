#!/usr/bin/env bash
set -euo pipefail
test -n "${SLURM_JOB_ID:-}" || { echo 'Use a Slurm compute allocation.' >&2;exit 2; }
project=$(cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$project"
threads=${BUILD_JOBS:-2}
test "$threads" -ge 1 && test "$threads" -le "${SLURM_CPUS_PER_TASK:-1}"
if [[ "${BUILD_RESUME_LINK_ONLY:-0}" == 1 ]];then
  sha256sum -c ../resume_sources_and_objects.sha256 > ../checked_resume_objects.log
else
  mkdir .native-u64
fi
build="$project/.native-u64"
source_dir="$project/Exam/NativeU64"
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
    make -B -C "${lib%/*}" -j "$jobs" "${common[@]}" "${lib##*/}"
    test -f "$lib"
done
/opt/ohpc/pub/cuda/13.0.2/bin/nvcc -O3 --use_fast_math --threads 1 -DXYZDBL -DGOTPM -DUSE_CUDA -DGFS_WEIGHT_DYNAMICS -DGFS_VARPRESS_W0 -DGFS_FIXEDPOLY_SIDES=1024 -gencode arch=compute_90,code=sm_90 "${inc[@]}" -c "$source_dir/exam_gpu.cu" -o "$build/exam_gpu.o"
/opt/ohpc/pub/cuda/13.0.2/bin/nvcc -O3 --use_fast_math --threads 1 -DXYZDBL -DGOTPM -DUSE_CUDA -DGFS_WEIGHT_DYNAMICS -DGFS_VARPRESS_W0 -DGFS_FIXEDPOLY_SIDES=1024 -gencode arch=compute_90,code=sm_90 "${inc[@]}" -c "$source_dir/gfs_unified_gpu.cu" -o "$build/gfs_unified_gpu.o"
mpiicx "${flags[@]}" "${inc[@]}" -c "$project/Exam/gfs_viscosity_implicit.c" -o "$build/gfs_viscosity_implicit.o"
for unit in exam exam_gpu_extract mpirks.exam2d exam3d_gfs gfs_circle_observer gfs_circle_dense kh util gfs_weight_checkpoint;do
    mpiicx "${flags[@]}" "${inc[@]}" -c "$source_dir/$unit.c" -o "$build/$unit.o"
done
fi
mpiicx "${flags[@]}" "${inc[@]}" -I"$project/FFT" -I"$project/Params" -I"$project/SPH" -I"$project/ParallelIO" -I"$project/Exam/RT" -I"$project/Exam/Kepler" -I"$project/Exam/MkGlass2D" -I"$project/Exam/Cylinder" -I"$project/Exam/Sedov2D" -c eunha2.c -o "$build/eunha2.o"
ar rcs "$build/libexam.a" "$build/exam.o" "$build/exam_gpu.o" "$build/gfs_unified_gpu.o" "$build/gfs_viscosity_implicit.o" "$build/exam_gpu_extract.o" "$build/mpirks.exam2d.o" "$build/exam3d_gfs.o" "$build/gfs_circle_observer.o" "$build/gfs_circle_dense.o"
ar rcs "$build/libkh.a" "$build/kh.o" "$build/util.o" "$build/gfs_weight_checkpoint.o"
# Preserve the qualified historical archive order; use only fresh local libraries.
libpaths=(-L"$build")
for lib in "${archives[@]}";do libpaths+=(-L"$project/${lib%/*}");done
libs=(-lheader -lcosmology -ltwolpt -lkjhrw -lrks -ltimerutil -lsph -lrecfast -ldenrw -lmpiaux -lparallelIO -lkh -lrt -lrtold -lkp -lglass2d -lcyl -lsedov2d -lexam -lost -lcolor -lmyvoro -lget_lambda -lnew_ssp -lextforce -lfft -lcamb)
mpiifx -nofor-main -g -O2 -qopenmp -DINTEL -DXYZDBL -DGOTPM "$build/eunha2.o" "${libpaths[@]}" "${libs[@]}" -L/home/kjhan/local/lib -lfftw3f_mpi -lfftw3f -lfftw3f_threads -L/opt/ohpc/pub/cuda/13.0.2/lib64 -lcudart -lcudadevrt -lstdc++ -o "$build/eunha2_native_u64.exe"
sha256sum "$build/eunha2_native_u64.exe" "$build/"*.o "$build/"*.a "${archives[@]}" > "$build/SHA256SUMS"
printf 'Full native U64 source build complete. Physical parity not yet qualified.\n' > "$build/COMPILE_COMPLETE"
