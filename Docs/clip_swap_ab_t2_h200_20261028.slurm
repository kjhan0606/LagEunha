#!/usr/bin/env bash
# Three paired, single-step timing and equivalence runs from the same native
# t=2 checkpoint. Baseline and candidate execute sequentially on one H200.

#SBATCH --job-name=KepClipSwapAB
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --chdir=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/diag_clip_swap_ab_t2_h200_20261028
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/diag_clip_swap_ab_t2_h200_20261028/slurm-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/diag_clip_swap_ab_t2_h200_20261028/slurm-%j.err

set -euo pipefail

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2

project=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026
q1="$project/q1"
run_dir="$q1/dteta_no_cap_from_t1p5_to65_20261028"
checkpoint="$run_dir/khout.016489.dat"
candidate="$PWD/candidate_binary.exe"
baseline="$run_dir/eunha2.exe"
python=/home/kjhan/miniconda3/bin/python3
mpi=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun

expected_baseline_sha=a6f16e8b1e143f70e8b8c32155008077bd84e813f7f4e4aff2f04b23dfd04572
expected_candidate_sha=6718551061d70e24e09f0fe000a35b5cf1ac8c1a15b579f1e97e6d907b42baa0
actual_baseline_sha=$(sha256sum "$baseline" | awk '{print $1}')
actual_candidate_sha=$(sha256sum "$candidate" | awk '{print $1}')
[[ "$actual_baseline_sha" == "$expected_baseline_sha" ]] || {
  echo "baseline executable hash changed: $actual_baseline_sha" >&2; exit 2;
}
[[ "$actual_candidate_sha" == "$expected_candidate_sha" ]] || {
  echo "candidate executable hash mismatch: $actual_candidate_sha" >&2; exit 2;
}

snapshot_times=$("$python" - "$checkpoint" <<'PY'
import os
import struct
import sys

path = sys.argv[1]
size = os.path.getsize(path)
with open(path, "rb") as stream:
    magic, n, stride, stress = struct.unpack("<4i", stream.read(16))
    if magic != -1 or n != 147456 or stride != 576 or stress != 1:
        raise SystemExit("unexpected native checkpoint header")
    if size != 16 + n * stride + 52:
        raise SystemExit("incomplete native checkpoint")
    stream.seek(16 + n * stride)
    metadata = struct.unpack("<ddfiiddff", stream.read(52))
    time = metadata[0]
    if abs(time - 2.0) > 1.0e-3:
        raise SystemExit(f"expected t=2 checkpoint, got t={time}")
print(f"{time:.15f} {time + 1.0e-5:.15f}")
PY
)
read -r start_time stop_time <<< "$snapshot_times"
checkpoint_sha=$(sha256sum "$checkpoint" | awk '{print $1}')
restart_step=16489

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

printf 'checkpoint=%s\ncheckpoint_sha256=%s\nstart_time=%s\nstop_time=%s\nrestart_step=%s\n' \
  "$checkpoint" "$checkpoint_sha" "$start_time" "$stop_time" "$restart_step" > run-manifest.txt
printf 'baseline_sha256=%s\ncandidate_sha256=%s\n' \
  "$actual_baseline_sha" "$actual_candidate_sha" >> run-manifest.txt

run_one() {
  local name=$1 executable=$2 directory="$PWD/$1"
  mkdir "$directory"
  ln -s "$run_dir/params.dat" "$directory/params.dat"
  ln -s "$q1/khout.000000.dat" "$directory/khout.000000.dat"
  ln -s "$checkpoint" "$directory/khout.016489.dat"
  ln -s "$executable" "$directory/eunha2.exe"
  printf 'variant=%s start_time=%s stop_time=%s checkpoint_sha256=%s\n' \
    "$name" "$start_time" "$stop_time" "$checkpoint_sha" > "$directory/run-input.txt"
  local began=$SECONDS rc
  set +e
  (cd "$directory" && "$mpi" -np 1 ./eunha2.exe params.dat "$restart_step" > log 2>&1)
  rc=$?
  set -e
  printf 'EXIT:%s\n' "$rc" >> "$directory/log"
  printf '%s wall_seconds=%s exit=%s\n' "$name" "$((SECONDS-began))" "$rc" | tee -a run-manifest.txt
  [[ $rc -eq 0 ]] || return "$rc"
  [[ -s "$directory/khout.016490.dat" ]] || {
    echo "$name did not produce khout.016490.dat" >&2; return 3;
  }
}

# Interleave order to reduce first-run/clock-state bias in the paired timing.
run_one baseline1 "$baseline"
run_one candidate1 "$candidate"
run_one candidate2 "$candidate"
run_one baseline2 "$baseline"
run_one baseline3 "$baseline"
run_one candidate3 "$candidate"

for pair in 1 2 3; do
  if cmp -s "baseline${pair}/khout.016490.dat" "candidate${pair}/khout.016490.dat"; then
    echo "PAIR${pair}_OUTPUT_BYTEWISE_IDENTICAL=1" | tee -a run-manifest.txt
  else
    echo "PAIR${pair}_OUTPUT_BYTEWISE_IDENTICAL=0" | tee -a run-manifest.txt
    sha256sum "baseline${pair}/khout.016490.dat" "candidate${pair}/khout.016490.dat" | tee -a run-manifest.txt
    exit 4
  fi
done

for prefix in baseline candidate; do
  for pair in 2 3; do
    previous=$((pair-1))
    cmp -s "${prefix}${previous}/khout.016490.dat" "${prefix}${pair}/khout.016490.dat" || {
      echo "${prefix} output is nondeterministic across repetitions" >&2; exit 5;
    }
  done
done

for name in baseline1 candidate1 candidate2 baseline2 baseline3 candidate3; do
  echo "===== $name ====="
  rg '^\[GPU tess device profile\]|^\[GPU tess\]|^\[GPU tess profile\]|^\[RK4E\]|^Time=' "$name/log" || true
done
echo "AB_RESULT=BYTEWISE_EQUIVALENT"
