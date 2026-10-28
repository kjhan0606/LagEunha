#!/usr/bin/env bash
#SBATCH --job-name=LagQ2ClipEarly
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A10:1
#SBATCH --mem=32G
#SBATCH --time=00:30:00
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
{ hostname; nvidia-smi -L; printf 'CUDA_VISIBLE_DEVICES=%s\n' "${CUDA_VISIBLE_DEVICES:-unset}"; } > device.log
sha256sum -c inputs.sha256
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0.06
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP=0.835 EUNHA_DUMP_DT=0.825 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_KEPLER_HALO_CELLS=32
export GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
export GFS_GPU_TESS_PROFILE=1
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS LAGMFM_GLASS_FILE EUNHA_KEPLER_INITIAL_AREAS_FILE
unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
TIMEFORMAT='[RUN_TIME] wall=%R'
set +e
time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat 426 > log 2>&1
run_rc=$?
set -e
printf 'EXIT:%s\n' "$run_rc" >> log
if [ "$run_rc" -ne 0 ]; then exit "$run_rc"; fi
python=/home/kjhan/miniconda3/bin/python3
endpoint=$($python - <<'PY'
from pathlib import Path
import sys
sys.path.insert(0, '/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs')
from prepare_kepler_sparse_area_restart_20261027 import snapshot_time
matches=[p for p in Path('.').glob('khout.*.dat') if abs(snapshot_time(p,55931)-0.835)<1e-10]
if len(matches)!=1:
    raise SystemExit(f'Expected one complete t=0.835 endpoint, found {matches}')
print(matches[0])
PY
)
initial=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_wide_area_t001/khout.000000.dat
$python audit_kepler_sparse_state_consistency_20261027.py "$initial" "$endpoint" --expected-count 55931 --time 0.835 --output audit_state_t0835.json > audit_state.log
$python audit_kepler_native_radial_20261026.py "$initial" "$endpoint" --mode baseline --output audit_radial_t0835.json > audit_radial.log
$python audit_kepler_sparse_native_ballistic_20261027.py "$initial" "$endpoint" --expected-count 55931 --time 0.835 --output audit_native_ballistic_t0835.json > audit_geometry.log
$python audit_kepler_seam_provenance_20261027.py "$initial" "$endpoint" --time 0.835 --output audit_seam_t0835.json > audit_seam.log
