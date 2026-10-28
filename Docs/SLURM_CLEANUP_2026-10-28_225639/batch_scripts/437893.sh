#!/usr/bin/env bash
#SBATCH --job-name=Q2MPI4t2231
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:H200:4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --nice=10000
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
PARENT=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825/to2p221_h28
INITIAL=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/khout.000000.dat
PYTHON=/home/kjhan/miniconda3/bin/python3
cd "$PWD"
sha256sum -c scripts.sha256
restart_step=$($PYTHON -c 'import json; print(json.load(open("restart_plan.json"))["step"])')
sha256sum eunha2.exe params.dat "$(printf 'khout.%06d.dat' "$restart_step")" > inputs.sha256
sha256sum -c inputs.sha256
mapfile -t visible_gpus < <(nvidia-smi -L)
printf 'visible GPU count=%s; MPI ranks=4; OMP=2\n' "${#visible_gpus[@]}"
printf '%s\n' "${visible_gpus[@]}"
if [[ ${#visible_gpus[@]} -ne 4 ]]; then
  echo "Expected 4 visible GPUs, found ${#visible_gpus[@]}" >&2
  exit 2
fi
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
printf 'MPI scaling: ranks=4, interval=2.221->2.231, restart_step=%s\n' "$restart_step"
TIMEFORMAT='[RUN_TIME] wall=%R'
set +e
time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 4 ./eunha2.exe params.dat "$restart_step" > log 2>&1
run_rc=$?
set -e
printf 'EXIT:%s\n' "$run_rc" >> log
if [ "$run_rc" -ne 0 ]; then exit "$run_rc"; fi
"$PYTHON" "$PARENT/audit_kepler_polar_long_gate_20261027.py" "$INITIAL" "$PWD" --time 2.231 --output audit_t2p231.json > audit.log
endpoint=$($PYTHON -c 'import json,sys; print(json.load(open(sys.argv[1]))["endpoint"])' audit_t2p231.json)
"$PYTHON" "$PARENT/audit_kepler_native_radial_20261026.py" "$INITIAL" "$endpoint" --mode baseline --output audit_radial_t2p231.json > audit_radial.log
