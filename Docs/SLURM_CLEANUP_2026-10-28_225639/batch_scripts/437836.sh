#!/usr/bin/env bash
#SBATCH --job-name=Q2OMPBench
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --nice=10000
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -euo pipefail

case "${OMP_BENCH:-}" in
  2|8) ;;
  *) echo 'OMP_BENCH must be 2 or 8' >&2; exit 2 ;;
esac

BENCH_ROOT=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825/to2p221_h28/omp_thread_pair_t2221_to2231_20261027
PARENT=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825/to2p221_h28
INITIAL=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/khout.000000.dat
DEST="$BENCH_ROOT/omp${OMP_BENCH}"
PYTHON=/home/kjhan/miniconda3/bin/python3
MPIEXEC=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun

cd "$BENCH_ROOT"
sha256sum -c scripts.sha256
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS="$OMP_BENCH" OMP_DYNAMIC=FALSE OMP_STACKSIZE=256M OMP_DISPLAY_ENV=VERBOSE
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm

"$PYTHON" "$PARENT/prepare_kepler_polar_long_restart_20261027.py" \
  "$PARENT" "$DEST" "$INITIAL" --receipt audit_t2p221.json \
  --expected-time 2.221 --target-time 2.231 > "$DEST/prepare.log"
restart_step=$("$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))["step"])' "$DEST/restart_plan.json")
cd "$DEST"
sha256sum eunha2.exe params.dat "$(printf 'khout.%06d.dat' "$restart_step")" > inputs.sha256
sha256sum -c inputs.sha256

export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0.06
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP=2.231 EUNHA_DUMP_DT=0.5 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_KEPLER_HALO_CELLS=28
export GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS LAGMFM_GLASS_FILE
unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
unset GFS_RK4_DIRECT_U GFS_THERMAL64 EUNHA_KEPLER_INITIAL_AREAS_FILE

printf 'OMP benchmark: threads=%s, interval=2.221->2.231, restart_step=%s\n' "$OMP_BENCH" "$restart_step"
nvidia-smi -L
TIMEFORMAT='[RUN_TIME] wall=%R'
set +e
time "$MPIEXEC" -np 1 ./eunha2.exe params.dat "$restart_step" > log 2>&1
run_rc=$?
set -e
printf 'EXIT:%s\n' "$run_rc" >> log
if [ "$run_rc" -ne 0 ]; then exit "$run_rc"; fi

"$PYTHON" "$PARENT/audit_kepler_polar_long_gate_20261027.py" "$INITIAL" "$DEST" \
  --time 2.231 --output audit_t2p231.json > audit.log
endpoint=$("$PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1]))["endpoint"])' audit_t2p231.json)
"$PYTHON" "$PARENT/audit_kepler_native_radial_20261026.py" "$INITIAL" "$endpoint" \
  --mode baseline --output audit_radial_t2p231.json > audit_radial.log
