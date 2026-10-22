#!/usr/bin/env bash
#SBATCH --job-name=material-full-ref
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=00:30:00
set -euo pipefail
dest=$1
base=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
src=/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src
source "$dest/run_env.sh"
export OMP_NUM_THREADS=4 OMP_PROC_BIND=false I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 I_MPI_CXX=g++
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
export GFS_MATERIAL_PREPARE_REFERENCE=1 GFS_MATERIAL_BOUNDARY=material_constant_Pext
export GFS_GPU_RAW_FACES=1 GFS_GPU_FINAL_VOLUME=1 GFS_GPU_FACE_DUMP_CALL=0 GFS_GPU_TESS_DUMP=0
export HYDRO_TSTOP=10 EUNHA_DUMP_DT=100 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
unset GFS_DG_PROJECT_ONLY GFS_MOMENT_INIT_ONLY GFS_DG_ENRICHED GFS_SSPRK3
mkdir "$dest/run"
cp "$base/core_refinement/fine_integrated/khout.000001.dat" "$dest/run/khout.000001.dat"
cp "$base/core_refinement/fine_integrated/params.dat" "$dest/run/params.dat"
sha256sum -c "$dest/source.sha256" > "$dest/prelaunch_source_integrity.txt"
sha256sum -c "$dest/retained_input.sha256" > "$dest/prelaunch_input_integrity.txt"
printf 'BUILDING_ISOLATED_NATIVE\n' > "$dest/status.txt"
includes=(-I"$dest/Exam" -I"$src" -I"$src/Exam" -I"$src/Voro" -I/home/kjhan/local/include -I"$src/RkS" -I"$src/Cosmos" -I"$src/Timer" -I"$src/MpiAux" -I"$src/Colorize" -I"$src/OST")
set +e
if test -f "$dest/retained_build_objects.sha256"; then
    sha256sum -c "$dest/retained_build_objects.sha256" > "$dest/retained_build_integrity.txt"
    host=0;cuda=0
else
mpiicx -g -O2 -fno-strict-aliasing -qopenmp -DINTEL -DXYZDBL -DGOTPM -DUSE_CUDA "${includes[@]}" -MMD -MF "$dest/Exam/exam.d" -c "$dest/Exam/exam.c" -o "$dest/Exam/exam.o" > "$dest/build_exam.log" 2>&1
host=$?
nvcc -O3 --use_fast_math -DXYZDBL -DGOTPM -DUSE_CUDA -gencode arch=compute_90,code=sm_90 -gencode arch=compute_90,code=compute_90 "${includes[@]}" -MMD -MF "$dest/Exam/exam_gpu.d" -c "$dest/Exam/exam_gpu.cu" -o "$dest/Exam/exam_gpu.o" > "$dest/build_cuda.log" 2>&1
cuda=$?
fi
mpicxx -O2 -std=c++17 -fno-fast-math -ffp-contract=off -DXYZDBL -DGOTPM -DUSE_CUDA "${includes[@]}" -I/opt/ohpc/pub/cuda/13.0.2/include -MMD -MF "$dest/Exam/material_full_reference.d" -c "$dest/Exam/material_full_reference.cpp" -o "$dest/Exam/material_full_reference.o" > "$dest/build_material.log" 2>&1
material=$?
set -e
printf '{"native_host":%d,"native_CUDA":%d,"material_host":%d}\n' "$host" "$cuda" "$material" > "$dest/build_exit.json"
if test "$host" != 0 || test "$cuda" != 0 || test "$material" != 0; then printf 'BUILD_FAILED\n' > "$dest/status.txt"; exit 1; fi
cp "$src/Exam/libexam.a" "$dest/Exam/libexam.a"
ar rcs "$dest/Exam/libexam.a" "$dest/Exam/exam.o" "$dest/Exam/exam_gpu.o" "$dest/Exam/material_full_reference.o"
# Reuse the proven original link recipe, with isolated exam archive/output.
python3 - "$dest" "$src" "$base/build_link.log" <<'PY'
import pathlib,shlex,subprocess,sys
d,src,log=map(pathlib.Path,sys.argv[1:]);args=shlex.split(log.read_text().splitlines()[0]);i=args.index('-o');args[i+1]=str(d/'eunha2_material_reference.exe');args=[('-L'+str(d/'Exam')) if a=='-L./Exam' else a for a in args];args+=['-lquadmath','-lgcc_s']
(d/'link_command.txt').write_text(shlex.join(args)+'\n')
with (d/'build_link.log').open('w') as out:subprocess.run(args,cwd=src,stdout=out,stderr=subprocess.STDOUT,check=True)
PY
sha256sum "$dest/eunha2_material_reference.exe" > "$dest/binary.sha256"
ldd "$dest/eunha2_material_reference.exe" > "$dest/ldd.txt"
if grep -q 'not found' "$dest/ldd.txt"; then printf 'DYNAMIC_LIBRARY_MISSING\n' > "$dest/status.txt"; exit 2; fi
nvidia-smi --query-gpu=index,uuid,name,memory.total,memory.used --format=csv > "$dest/device.txt"
printf 'RUNNING_ORIGINAL_NATIVE_FULL_REFERENCE\n' > "$dest/status.txt"
cd "$dest/run"
env | sort > "$dest/effective_environment.txt"
set +e
mpirun -np 4 "$dest/eunha2_material_reference.exe" params.dat 1 > "$dest/native.log" 2>&1
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/retained_input.sha256" > "$dest/final_input_integrity.txt"
if test "$result" = 0; then printf 'FULL_REFERENCE_PREPARED_NO_EVOLUTION\n' > "$dest/status.txt"; else printf 'FULL_REFERENCE_REJECTED_EXIT_%s\n' "$result" > "$dest/status.txt"; fi
exit "$result"
