#!/usr/bin/env bash
#SBATCH --job-name=LagQ2W64FastBuild
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:0
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --nice=10000
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
sha256sum -c inputs.sha256
run_dir=$PWD
build_dir=${SLURM_TMPDIR:-/tmp}/lag_q2_fast_capacity64_${SLURM_JOB_ID}
mkdir -p "$build_dir"
tar -xf source.tar -C "$build_dir"
python3 patch_guard.py "$build_dir" > patch_guard_receipt.json
python3 patch_iterations.py "$build_dir" > patch_iterations_receipt.json
python3 patch_zero.py "$build_dir" > patch_zero_receipt.json
python3 patch_capacity.py "$build_dir" > patch_capacity_receipt.json
for patch_name in gap.patch accel.patch diagnostic.patch thermal.patch; do
    patch --batch --fuzz=0 --dry-run -d "$build_dir" -p1 < "$run_dir/$patch_name"
    patch --batch --fuzz=0 -d "$build_dir" -p1 < "$run_dir/$patch_name"
done > patch_performance.log
sha256sum "$build_dir/Exam/exam.c" > exam_source.sha256
g++ -std=c++14 -O2 -Wall -Wextra -Werror -I "$build_dir/Exam" test_projection.cpp -o "$build_dir/test_projection"
"$build_dir/test_projection" > test_projection.log
rm -f "$build_dir/Exam/exam.o" "$build_dir/Exam/exam_gpu.o" \
      "$build_dir/Exam/exam_gpu_extract.o" "$build_dir/Exam/gfs_energy_neutral_probe.o" \
      "$build_dir/Exam/libexam.a" "$build_dir/eunha2.exe"
cd "$build_dir"
TIMEFORMAT='[BUILD_TIME] wall=%R'
time make -j2 USE_CUDA=1 GFS_VARPRESS_W0=1 > "$run_dir/build.log" 2>&1
test -s eunha2.exe
cp eunha2.exe "$run_dir/eunha2_fast_capacity64_v6.exe"
sha256sum "$run_dir/eunha2_fast_capacity64_v6.exe" > "$run_dir/binary.sha256"
