#!/usr/bin/env bash
set -euo pipefail
dest=$(cd -- "$(dirname -- "$0")" && pwd)
src=/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src
test "$(hostname)" = LagEunha
source "$dest/run_env.sh"
export I_MPI_CXX=g++ OMP_PROC_BIND=false
sha256sum -c "$dest/source.sha256" > "$dest/prebuild_source_integrity.txt"
sha256sum -c "$dest/retained_build_objects.sha256" > "$dest/retained_build_integrity.txt"
printf 'BUILDING_NATIVE_PAIRED_TRANSITION_NO_GPU_EXECUTION\n' > "$dest/status.txt"
trap 'printf "BUILD_REJECTED_NO_PHYSICAL_EXECUTION\n" > "$dest/status.txt"' ERR
includes=(-I"$dest/Exam" -I"$src" -I"$src/Exam" -I"$src/Voro" -I/home/kjhan/local/include -I"$src/RkS" -I"$src/Cosmos" -I"$src/Timer" -I"$src/MpiAux" -I"$src/Colorize" -I"$src/OST")
nvcc -O2 -std=c++17 --fmad=false -DXYZDBL -DGOTPM -DUSE_CUDA -gencode arch=compute_89,code=sm_89 "${includes[@]}" -c "$dest/Exam/material_full_initial_cuda.cu" -o "$dest/Exam/material_full_initial_cuda.o" > "$dest/build_cuda.log" 2>&1
mpicxx -O2 -std=c++17 -fopenmp -fno-fast-math -ffp-contract=off -DXYZDBL -DGOTPM -DUSE_CUDA "${includes[@]}" -I/opt/ohpc/pub/cuda/13.0.2/include -c "$dest/Exam/material_full_initial_driver.cpp" -o "$dest/Exam/material_full_initial_driver.o" > "$dest/build_host.log" 2>&1
ar rcs "$dest/Exam/libexam.a" "$dest/Exam/material_full_initial_driver.o" "$dest/Exam/material_full_initial_cuda.o"
python3 - "$dest" "$src" <<'PY'
import pathlib,shlex,subprocess,sys
d,src=map(pathlib.Path,sys.argv[1:]); args=shlex.split((d/'donor_link_command.txt').read_text());args[args.index('-o')+1]=str(d/'eunha2_positive_transition.exe')
args=['-L'+str(d/'Exam') if a.startswith('-L') and a.endswith('/material_full_GPU_mass_v1_h200_mpi4/Exam') else a for a in args]
assert '-L'+str(d/'Exam') in args
(d/'link_command.txt').write_text(shlex.join(args)+'\n')
with (d/'build_link.log').open('x') as log:subprocess.run(args,cwd=src,stdout=log,stderr=subprocess.STDOUT,check=True)
PY
sha256sum "$dest/eunha2_positive_transition.exe" "$dest/Exam/libexam.a" > "$dest/binary.sha256"
ldd "$dest/eunha2_positive_transition.exe" > "$dest/ldd.txt"
if rg -q 'not found' "$dest/ldd.txt"; then exit 2; fi
sha256sum -c "$dest/source.sha256" > "$dest/postbuild_source_integrity.txt"
printf 'BUILD_PASS_NO_PHYSICAL_EXECUTION\n' > "$dest/status.txt"
