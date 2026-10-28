#!/usr/bin/env bash
#SBATCH --job-name=LagQ2H28To30
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=36:00:00
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
predecessor_time=28; target_time=30; predecessor_receipt=audit_t28.json; tag=t30
sha256sum -c scripts.sha256
parent=$(dirname "$PWD")
initial=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/khout.000000.dat
python=/home/kjhan/miniconda3/bin/python3
printf '2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43  %s\n' "$parent/params.dat" | sha256sum -c -
"$python" prepare_kepler_polar_long_restart_20261027.py "$parent" "$PWD" "$initial" \
  --receipt "$predecessor_receipt" --expected-time "$predecessor_time" \
  --target-time "$target_time" > prepare.log
restart_step=$($python -c 'import json; print(json.load(open("restart_plan.json"))["step"])')
sha256sum eunha2.exe params.dat "$(printf 'khout.%06d.dat' "$restart_step")" > inputs.sha256
sha256sum -c inputs.sha256
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0.06
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP="$target_time" EUNHA_DUMP_DT=0.5 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_KEPLER_HALO_CELLS=28
export GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS LAGMFM_GLASS_FILE
unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
unset GFS_RK4_DIRECT_U GFS_THERMAL64 EUNHA_KEPLER_INITIAL_AREAS_FILE
TIMEFORMAT='[RUN_TIME] wall=%R'
set +e
time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat "$restart_step" > log 2>&1
run_rc=$?
set -e
printf 'EXIT:%s\n' "$run_rc" >> log
if [ "$run_rc" -ne 0 ]; then exit "$run_rc"; fi
"$python" audit_kepler_polar_long_gate_20261027.py "$initial" "$PWD" \
  --time "$target_time" --output "audit_${tag}.json" > audit.log
endpoint=$($python -c 'import json,sys; print(json.load(open(sys.argv[1]))["endpoint"])' "audit_${tag}.json")
"$python" audit_kepler_native_radial_20261026.py "$initial" "$endpoint" \
  --mode baseline --output "audit_radial_${tag}.json" > audit_radial.log
