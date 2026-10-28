#!/usr/bin/env bash
#SBATCH --job-name=KepAB
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=0-12:00:00
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=8 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
PY=/home/kjhan/miniconda3/bin/python3
PARENT=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/dteta_no_cap_from_t1p5_to65_20261028
ROOT=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1
if [[ "${INTEGRATOR:-}" == 1 ]]; then
  grep -q '^define Hydro time-stepping way = 1 # RK4 control$' params.dat
elif [[ "${INTEGRATOR:-}" == 3 ]]; then
  grep -q '^define Hydro time-stepping way = 3 # KDK blend$' params.dat
else
  echo "INTEGRATOR must be 1 (RK4) or 3 (KDK blend)" >&2
  exit 2
fi
sha256sum -c static_inputs.sha256
for f in log restart.sha256 audit_t6p5_restart_gate.json audit_t7_gate.json audit_t7_radial.json; do
  if [[ -e "$f" ]]; then echo "Refusing to overwrite $f" >&2; exit 2; fi
done
"$PY" ./audit_kepler_polar_long_gate_20261027.py "$ROOT/khout.000000.dat" "$PARENT" --time 6.5 --output audit_t6p5_restart_gate.json > audit_t6p5_restart_gate.log
ENDPOINT=$("$PY" -c 'import json; print(json.load(open("audit_t6p5_restart_gate.json"))["endpoint"])')
STEP=$(basename "$ENDPOINT" | sed -E 's/^khout\.([0-9]+)\.dat$/\1/')
if [[ ! "$STEP" =~ ^[0-9]+$ ]]; then echo "invalid t=6.5 restart endpoint: $ENDPOINT" >&2; exit 2; fi
cp "$ENDPOINT" "khout.${STEP}.dat"
sha256sum "khout.${STEP}.dat" > restart.sha256
sha256sum -c restart.sha256
printf '%s\n' "$STEP" > restart_step.txt
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP=7.0 EUNHA_DUMP_DT=0.5 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_KEPLER_HALO_CELLS=28
export GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
TIMEFORMAT='[RUN_TIME] wall=%R'
set +e
time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat "$STEP" > log 2>&1
RC=$?
set -e
printf 'EXIT:%s\n' "$RC" >> log
if [[ "$RC" -ne 0 ]]; then exit "$RC"; fi
"$PY" ./audit_kepler_polar_long_gate_20261027.py "$ROOT/khout.000000.dat" "$PWD" --time 7.0 --output audit_t7_gate.json > audit_t7.log
ENDPOINT7=$("$PY" -c 'import json; print(json.load(open("audit_t7_gate.json"))["endpoint"])')
"$PY" ./audit_kepler_native_radial_20261026.py "$ROOT/khout.000000.dat" "$ENDPOINT7" --mode baseline --output audit_t7_radial.json > audit_t7_radial.log
sha256sum khout.*.dat audit_t7_gate.json audit_t7_radial.json > outputs.sha256
