#!/usr/bin/env bash
#SBATCH --job-name=LagPolarResume
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
parent_job_id=$(cat parent_job_id.txt)
[[ "$parent_job_id" =~ ^[0-9]+$ ]]
parent_state=$(sacct -X -j "$parent_job_id" --noheader --format=State | awk 'NF {print $1; exit}')
parent_state=${parent_state%%+}
printf '[PARENT] job=%s state=%s\n' "$parent_job_id" "$parent_state"
case "$parent_state" in
    TIMEOUT|PREEMPTED|NODE_FAIL|BOOT_FAIL) ;;
    *) echo 'Refusing automatic restart after a physical or user stop'; exit 21 ;;
esac
preflight=$(dirname "$PWD")
parent="$preflight/to0p825"
python=/home/kjhan/miniconda3/bin/python3
"$python" prepare_kepler_polar_timeout_resume_20261026.py \
    "$preflight" "$parent" "$PWD" --parent-state="$parent_state" > prepare.log
mode=$($python -c 'import json; print(json.load(open("resume_plan.json"))["mode"])')
restart_step=$($python -c 'import json; print(json.load(open("resume_plan.json"))["step"])')
checkpoint=$(printf 'khout.%06d.dat' "$restart_step")
sha256sum eunha2.exe params.dat "$checkpoint" > inputs.sha256
sha256sum -c inputs.sha256
if [[ "$mode" == resume ]]; then
    export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
    export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
    export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0.06
    export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
    export HYDRO_TSTOP=0.825 EUNHA_DUMP_DT=0.20625 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
    export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_KEPLER_HALO_CELLS=8
    export GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
    export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
    unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME LAGMFM_GLASS_FILE
    unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
    unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
    TIMEFORMAT='[RUN_TIME] wall=%R'
    set +e
    time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat "$restart_step" > log 2>&1
    run_rc=$?
    set -e
    printf 'EXIT:%s\n' "$run_rc" >> log
    if [ "$run_rc" -ne 0 ]; then exit "$run_rc"; fi
elif [[ "$mode" != audit_only ]]; then
    echo "Unexpected resume mode: $mode"
    exit 22
fi
"$python" audit_kepler_polar_t0825_20261026.py "$preflight/khout.000000.dat" "$PWD" > audit.log
