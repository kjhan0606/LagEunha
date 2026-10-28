#!/usr/bin/env bash
#SBATCH --job-name=LagSparseLongA10
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A10:1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
{ hostname; nvidia-smi -L; printf 'CUDA_VISIBLE_DEVICES=%s\n' "${CUDA_VISIBLE_DEVICES:-unset}"; } > device.log
if [ "$#" -ne 7 ]; then
  echo 'usage: runner parent count predecessor_time target_time dump_dt receipt consistency_receipt' >&2
  exit 2
fi
parent=$1
count=$2
predecessor_time=$3
target_time=$4
dump_dt=$5
receipt=$6
consistency_receipt=$7
python=/home/kjhan/miniconda3/bin/python3
"$python" - "$predecessor_time" "$target_time" <<'PY'
import math
import sys

try:
    predecessor, target = map(float, sys.argv[1:])
except ValueError as exc:
    raise SystemExit(f"invalid long-run time: {exc}")
if not (math.isfinite(predecessor) and math.isfinite(target)
        and 10 <= predecessor < target):
    raise SystemExit("long-run times must satisfy finite 10 <= predecessor < target")
PY
if [[ ! "$dump_dt" =~ ^[0-9]+$ ]] || [ "$dump_dt" -le 0 ]; then
  echo 'checkpoint interval must be a positive integer time' >&2
  exit 2
fi
sha256sum -c scripts.sha256
"$python" prepare_kepler_sparse_long_restart_20261027.py "$parent" "$PWD" \
  --receipt "$receipt" --consistency-receipt "$consistency_receipt" \
  --expected-count "$count" --expected-time "$predecessor_time" \
  --target-time "$target_time" > prepare.log
restart_step=$($python -c 'import json; print(json.load(open("restart_plan.json"))["step"])')
initial=$($python -c 'import json; print(json.load(open("restart_plan.json"))["initial"])')
sha256sum eunha2.exe params.dat "$(printf 'khout.%06d.dat' "$restart_step")" > inputs.sha256
sha256sum -c inputs.sha256
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0.06
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP="$target_time" EUNHA_DUMP_DT="$dump_dt" GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_KEPLER_HALO_CELLS=32
export GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
export GFS_GPU_TESS_PROFILE=1
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS LAGMFM_GLASS_FILE EUNHA_KEPLER_INITIAL_AREAS_FILE
unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
TIMEFORMAT='[RUN_TIME] wall=%R'
set +e
time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat "$restart_step" > log 2>&1
run_rc=$?
set -e
printf 'EXIT:%s\n' "$run_rc" >> log
if [ "$run_rc" -ne 0 ]; then exit "$run_rc"; fi
"$python" audit_kepler_sparse_area_t10_20261027.py "$initial" "$PWD" \
  --expected-count "$count" --time "$target_time" \
  --output "audit_t${target_time}.json" > audit.log
endpoint=$($python -c 'import json,sys; print(json.load(open(sys.argv[1]))["endpoint"])' "audit_t${target_time}.json")
"$python" audit_kepler_sparse_state_consistency_20261027.py "$initial" "$endpoint" \
  --expected-count "$count" --time "$target_time" \
  --output "audit_consistency_t${target_time}.json" > audit_consistency.log
"$python" audit_kepler_sparse_native_ballistic_20261027.py "$initial" "$endpoint" \
  --expected-count "$count" --time "$target_time" \
  --output "audit_ballistic_t${target_time}.json" > audit_ballistic.log
"$python" audit_kepler_seam_provenance_20261027.py "$initial" "$endpoint" \
  --time "$target_time" \
  --output "audit_seam_t${target_time}.json" > audit_seam.log
