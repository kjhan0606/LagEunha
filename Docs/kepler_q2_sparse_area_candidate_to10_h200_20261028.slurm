#!/usr/bin/env bash
# Resume the last complete sparse buffered-area Q2 state with the tested
# clip-swap CUDA candidate; stop at the exact t=10 validation gate.

#SBATCH --job-name=KepSparseClipTo10
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --chdir=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_buffered_area_t001/to0p825/to10_clip_swap_h200_20261028
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_buffered_area_t001/to0p825/to10_clip_swap_h200_20261028/slurm-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_buffered_area_t001/to0p825/to10_clip_swap_h200_20261028/slurm-%j.err

set -euo pipefail

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2

project=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026
source_run="$project/q2_profile_buffered_area_t001/to0p825/to10"
checkpoint="$source_run/khout.001318.dat"
initial="$project/q2_profile_buffered_area_t001/khout.000000.dat"
params="$source_run/params.dat"
candidate="$project/q1/diag_clip_swap_ab_t2_h200_20261028/candidate_binary.exe"
audit="$source_run/audit_kepler_sparse_area_t10_20261027.py"
python=/home/kjhan/miniconda3/bin/python3
mpi=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun

expected_checkpoint_sha=c3525f7486e143a65319e7d99debc2fff2240fd5a411354f4f6f80fe5c4e832a
expected_initial_sha=001dcd63e05937c9fabf0e12ba260dfa36369bed31baa762040936102523149b
expected_params_sha=2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43
expected_candidate_sha=6718551061d70e24e09f0fe000a35b5cf1ac8c1a15b579f1e97e6d907b42baa0

for output in params.dat eunha2.exe khout.000000.dat khout.001318.dat \
              log run-manifest.txt inputs.sha256 audit_t10.json audit_t10.log; do
  if [[ -e "$output" || -L "$output" ]]; then
    echo "Refusing to overwrite $output" >&2
    exit 2
  fi
done

for pair in \
  "$expected_checkpoint_sha $checkpoint" \
  "$expected_initial_sha $initial" \
  "$expected_params_sha $params" \
  "$expected_candidate_sha $candidate"; do
  read -r expected path <<< "$pair"
  actual=$(sha256sum "$path" | awk '{print $1}')
  [[ "$actual" == "$expected" ]] || {
    echo "Input hash mismatch: $path expected=$expected actual=$actual" >&2
    exit 2
  }
done

"$python" - "$checkpoint" <<'PY'
import os
import struct
import sys

path = sys.argv[1]
with open(path, "rb") as stream:
    header = struct.unpack("<4i", stream.read(16))
    stream.seek(-52, os.SEEK_END)
    time = struct.unpack("<d", stream.read(8))[0]
if header != (-1, 54312, 576, 1):
    raise SystemExit(f"Unexpected checkpoint header: {header}")
if os.path.getsize(path) != 16 + 54312 * 576 + 52:
    raise SystemExit("Checkpoint has an incomplete native record payload")
if abs(time - 2.4753977144699952) > 1e-12:
    raise SystemExit(f"Unexpected restart time: {time}")
print(f"checkpoint_valid=1 time={time:.15f} n=54312")
PY

cp --preserve=mode,timestamps "$params" params.dat
cp --preserve=mode,timestamps "$initial" khout.000000.dat
cp --preserve=mode,timestamps "$checkpoint" khout.001318.dat
cp --preserve=mode,timestamps "$candidate" eunha2.exe
sha256sum eunha2.exe params.dat khout.000000.dat khout.001318.dat > inputs.sha256
sha256sum -c inputs.sha256

cat > run-manifest.txt <<EOF
purpose=continue sparse buffered-area Q2 checkpoint with clip-swap tessellation candidate
source_run=$source_run
restart_checkpoint=$checkpoint
restart_checkpoint_sha256=$expected_checkpoint_sha
restart_step=1318
restart_time=2.4753977144699952
target_time=10.0
candidate_binary=$candidate
candidate_binary_sha256=$expected_candidate_sha
initial=$initial
initial_sha256=$expected_initial_sha
W=0
mass_mode=fixed
nominal_grid=384x384
actual_cell_count=54312
geometry=periodic-area sparse buffered disk
EOF

export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0.06
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP=10 EUNHA_DUMP_DT=0.5 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_KEPLER_HALO_CELLS=32
export GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=0.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS
unset LAGMFM_GLASS_FILE EUNHA_KEPLER_INITIAL_AREAS_FILE
unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION

TIMEFORMAT='[RUN_TIME] wall=%R'
set +e
time "$mpi" -np 1 ./eunha2.exe params.dat 1318 > log 2>&1
run_rc=$?
set -e
printf 'EXIT:%s\n' "$run_rc" >> log
printf 'exit=%s\n' "$run_rc" >> run-manifest.txt
[[ "$run_rc" -eq 0 ]] || exit "$run_rc"

if ! grep -q '\[GPU_CANDIDATE_STREAM\].*all_candidates_accounted_for=1' log; then
  echo "Candidate stream accounting marker is missing" >&2
  exit 3
fi
if grep -Eq '\[GPU_FINAL_VOLUME\].*invalid=[1-9]|application called MPI_Abort' log; then
  echo "Invalid GPU geometry or MPI abort marker found" >&2
  exit 4
fi

"$python" "$audit" "$PWD/khout.000000.dat" "$PWD" \
  --expected-count 54312 --output audit_t10.json > audit_t10.log
printf 'audit_t10=PASS\n' >> run-manifest.txt
