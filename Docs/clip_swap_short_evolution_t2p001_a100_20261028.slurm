#!/usr/bin/env bash
# Short matched evolution from the same t=2 checkpoint to expose whether the
# pointer-swap tessellator's near-degenerate face-count difference grows.

#SBATCH --job-name=KepClipSwapEvo
#SBATCH --partition=a100_pcie
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A100:1
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --chdir=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/diag_clip_swap_evolve_t2p001_a100_20261028
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/diag_clip_swap_evolve_t2p001_a100_20261028/slurm-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/diag_clip_swap_evolve_t2p001_a100_20261028/slurm-%j.err

set -euo pipefail

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2

project=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026
q1="$project/q1"
run_dir="$q1/dteta_no_cap_from_t1p5_to65_20261028"
checkpoint="$run_dir/khout.016489.dat"
candidate_binary="$q1/diag_clip_swap_ab_t2_a100_pcie_20261028/candidate_binary.exe"
python=/home/kjhan/miniconda3/bin/python3
mpi=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun
compare_script=/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs/compare_kepler_native_state_20261028.py

snapshot_times=$("$python" - "$checkpoint" <<'PY'
import os
import struct
import sys

path = sys.argv[1]
size = os.path.getsize(path)
with open(path, "rb") as stream:
    magic, n, stride, stress = struct.unpack("<4i", stream.read(16))
    if (magic, n, stride, stress) != (-1, 147456, 576, 1):
        raise SystemExit("unexpected native checkpoint header")
    if size != 16 + n * stride + 52:
        raise SystemExit("incomplete native checkpoint")
    stream.seek(16 + n * stride)
    metadata = struct.unpack("<ddfiiddff", stream.read(52))
    if abs(metadata[0] - 2.0) > 1.0e-3:
        raise SystemExit(f"expected t=2 checkpoint, got t={metadata[0]}")
print(f"{metadata[0]:.15f}")
PY
)
read -r start_time <<< "$snapshot_times"
restart_step=16489
stop_time=2.001

export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP="$stop_time" EUNHA_DUMP_DT=0.5 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_KEPLER_HALO_CELLS=28
export GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=0.1
export GFS_GPU_LOCAL_PERIODIC_SYNC=0 GFS_GPU_TESS_PROFILE=1
unset GFS_VARPRESS_GEOM_AUDIT GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS
unset LAGMFM_GLASS_FILE EUNHA_KEPLER_INITIAL_AREAS_FILE

printf 'checkpoint=%s\nstart_time=%s\nstop_time=%s\nrestart_step=%s\n' \
  "$checkpoint" "$start_time" "$stop_time" "$restart_step" > run-manifest.txt

run_one() {
  local name=$1 executable=$2 directory="$PWD/$1"
  mkdir "$directory"
  ln -s "$run_dir/params.dat" "$directory/params.dat"
  ln -s "$q1/khout.000000.dat" "$directory/khout.000000.dat"
  ln -s "$checkpoint" "$directory/khout.016489.dat"
  ln -s "$executable" "$directory/eunha2.exe"
  printf 'variant=%s start_time=%s stop_time=%s\n' \
    "$name" "$start_time" "$stop_time" > "$directory/run-input.txt"
  local began=$SECONDS rc
  set +e
  (cd "$directory" && "$mpi" -np 1 ./eunha2.exe params.dat "$restart_step" > log 2>&1)
  rc=$?
  set -e
  printf 'EXIT:%s\n' "$rc" >> "$directory/log"
  printf '%s wall_seconds=%s exit=%s\n' "$name" "$((SECONDS-began))" "$rc" | tee -a run-manifest.txt
  [[ $rc -eq 0 ]] || return "$rc"
  local step
  step=$(awk '$1 == "Time=" {value=$4} END {if (value != "") print value}' "$directory/log")
  [[ "$step" =~ ^[0-9]+$ ]] || { echo "$name has no final Time/step line" >&2; return 3; }
  local output="$directory/khout.$(printf '%06d' "$step").dat"
  [[ -s "$output" ]] || { echo "$name missing native endpoint $output" >&2; return 4; }
  printf '%s\n' "$output" > "$directory/final-output.txt"
  printf '%s step=%s output=%s\n' "$name" "$step" "$output" | tee -a run-manifest.txt
  rg '^\[GPU tess\]|^\[RK4E\]|^Time=' "$directory/log" || true
}

run_one baseline "$run_dir/eunha2.exe"
run_one candidate "$candidate_binary"

baseline_output=$(cat baseline/final-output.txt)
candidate_output=$(cat candidate/final-output.txt)
"$python" "$compare_script" "$baseline_output" "$candidate_output" \
  --output semantic_compare.json > semantic_compare.log
"$python" - semantic_compare.json <<'PY' | tee -a run-manifest.txt
import json
result = json.load(open("semantic_compare.json"))
print("primary_state_exact=" + str(result["primary_state_exact"]))
print("all_numeric_finite=" + str(result["all_numeric_finite"]))
print("numeric_difference_fields=" + ",".join(result["numeric_differences"]))
print("reference_t=" + str(result["metadata"]["reference"]["t"]))
print("candidate_t=" + str(result["metadata"]["candidate"]["t"]))
if not result["ids_equal"] or not result["all_numeric_finite"]:
    raise SystemExit("native state comparison failed structural/finiteness checks")
PY
