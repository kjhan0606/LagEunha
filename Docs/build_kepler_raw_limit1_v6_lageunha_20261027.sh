#!/usr/bin/env bash
# Incremental, isolated build on the user-authorized lageunha CPU server.
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
docs=/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs
root=/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027
parent_tree=/tmp/lag_q2_energy_neutral_inplace_v5_manual_001
build_tree=/tmp/lag_q2_energy_neutral_raw_limit1_v6_manual_001
output=$root/gpu_energy_neutral_native_build_v6_raw_limit1
test -d "$parent_tree"
test ! -e "$build_tree"
test ! -e "$output"
mkdir "$output"
cp "$docs/patch_kepler_raw_limit1_specialization_v6_20261027.py" "$output/patch_specialize.py"
sha256sum "$output/patch_specialize.py" > "$output/patch_specialize.sha256"
mv "$parent_tree" "$build_tree"
python3 "$output/patch_specialize.py" "$build_tree" > "$output/patch_specialize_receipt.json"
cd "$build_tree"
TIMEFORMAT='[BUILD_TIME] wall=%R'
time make -j2 USE_CUDA=1 GFS_VARPRESS_W0=1 > "$output/build.log" 2>&1
test -s eunha2.exe
cp eunha2.exe "$output/eunha2_energy_neutral_raw_limit1_v6.exe"
sha256sum "$output/eunha2_energy_neutral_raw_limit1_v6.exe" > "$output/binary.sha256"
tail -n 14 "$output/build.log"
