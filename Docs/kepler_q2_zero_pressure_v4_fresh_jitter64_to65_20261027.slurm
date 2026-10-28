#!/usr/bin/env bash
#SBATCH --job-name=LagQ2W64Fresh65v4
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A10:1
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --nice=200
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

pilot=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_zero_pressure_v4_t001_a10_v1')
run=Path.cwd()
receipt=json.loads((pilot/'audit_t001.json').read_text())
pair=json.loads((pilot/'pair_audit.json').read_text())
if not all(receipt['checks'].values()) or abs(receipt['time']-.001)>1e-12:
    raise ValueError('Fresh weighted pilot did not pass its native gate')
if not all(pair['checks'].values()) or pair['v4_endpoint_sha256']!=receipt['endpoint_sha256']:
    raise ValueError('v4 native regression did not pass on this pilot endpoint')
if (receipt['reference_sha256']!='8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b'
        or receipt['binary_sha256']!='3120aea6dad944b144e04d4597cb3f644abc6c803ff2318b5862b6bf8eec0adc'):
    raise ValueError('Fresh pilot lineage differs from staged inputs')
source=pilot/'khout.000001.dat'
if digest(source)!=receipt['endpoint_sha256'] or abs(snapshot_time(source,55351)-.001)>1e-10:
    raise ValueError('Fresh pilot checkpoint changed or is incomplete')
target=run/'khout.000001.dat'
if target.exists():
    raise FileExistsError(target)
shutil.copy2(source,target)
if digest(target)!=receipt['endpoint_sha256']:
    raise ValueError('Copied restart differs from accepted pilot')
with (run/'restart_origin.json').open('x') as output:
    json.dump({'pilot_job':437452,'pilot_receipt':str(pilot/'audit_t001.json'),
               'pair_receipt':str(pilot/'pair_audit.json'),
               'restart_sha256':receipt['endpoint_sha256'],'time':.001},output,indent=2)
    output.write('\n')
PY_RESTART
export GFS_GPU_VOLW2=1 GFS_GPU_VOLW2_ITERS=64
export GFS_GPU_VOLW2_MASK="$PWD/disk_active_mask.bin"
export GFS_KEPLER_FIXED_POLYGON=1 GFS_KEPLER_HALO_CELLS=1
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0.06
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP=6.5 EUNHA_DUMP_DT=0.825 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
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
    initial.dat khout.000001.dat . --parent-time .001 --time 6.5 \
    --expected-binary-sha256 3120aea6dad944b144e04d4597cb3f644abc6c803ff2318b5862b6bf8eec0adc \
    --output audit_t65.json > audit.log 2>&1
