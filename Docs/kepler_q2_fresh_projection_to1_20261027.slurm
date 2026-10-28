#!/usr/bin/env bash
#SBATCH --job-name=LagQ2Fresh__MODE__T1
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A10:1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --nice=20000
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
sha256sum -c inputs.sha256
/home/kjhan/miniconda3/bin/python3 - <<'PY_RESTART'
import json
from pathlib import Path
import shutil
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time

parent=Path('__PARENT_DIR__')
run=Path.cwd()
receipt=json.loads((parent/'audit_t001.json').read_text())
if (receipt['time']!=.001 or len(receipt['checks'])!=18
        or not all(receipt['checks'].values())
        or receipt['binary_sha256']!='__BINARY_SHA__'
        or receipt['reference_sha256']!='__INITIAL_SHA__'):
    raise ValueError('Fresh one-step parent failed provenance or physics gate')
source=parent/'khout.000001.dat'
if digest(source)!=receipt['endpoint_sha256'] or abs(snapshot_time(source,55351)-.001)>1e-10:
    raise ValueError('Fresh one-step parent checkpoint changed')
target=run/'khout.000001.dat'
if target.exists():
    raise FileExistsError(target)
shutil.copy2(source,target)
if digest(target)!=receipt['endpoint_sha256']:
    raise ValueError('Copied fresh parent differs')
with (run/'restart_origin.json').open('x') as output:
    json.dump({'parent_receipt':str(parent/'audit_t001.json'),
               'parent_sha256':receipt['endpoint_sha256'],'time':.001},output,indent=2)
    output.write('\n')
PY_RESTART
export GFS_GPU_VOLW2=1 GFS_GPU_VOLW2_ITERS=64 GFS_GPU_VOLW2_PROJECTION=__PROJECTION__
export GFS_GPU_VOLW2_MASK="$PWD/disk_active_mask.bin"
export GFS_KEPLER_FIXED_POLYGON=1 GFS_KEPLER_HALO_CELLS=1
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0.06
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP=1.0 EUNHA_DUMP_DT=0.1 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS LAGMFM_GLASS_FILE EUNHA_KEPLER_INITIAL_AREAS_FILE
unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
unset GFS_GPU_WDOT_PROBE_MASK GFS_GPU_FACE_DUMP_CALL GFS_GPU_FACE_DUMP_IDS GFS_GPU_TESS_DUMP
nvidia-smi -L
TIMEFORMAT='[RUN_TIME] wall=%R'
time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat 1 > log 2>&1
/home/kjhan/miniconda3/bin/python3 audit_kepler_q2_weighted_equilibrium_20261027.py \
    initial.dat khout.000001.dat . --parent-time .001 --time 1.0 \
    --expected-binary-sha256 __BINARY_SHA__ \
    --output audit_t1.json > audit.log 2>&1
