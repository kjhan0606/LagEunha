#!/usr/bin/env bash
#SBATCH --job-name=LagQ2ThermPair
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --nice=100
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=8 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
root=$PWD
sha256sum -c inputs.sha256
/home/kjhan/miniconda3/bin/python3 - <<'PY_PARENT'
import json
from pathlib import Path
import shutil
from prepare_kepler_sparse_area_restart_20261027 import digest, snapshot_time

root=str(Path.cwd())
parent=Path('/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_capacity64_v5_fresh_to30_omp8_h200_v1')
receipt_path=parent/'audit_t30.json'
receipt=json.loads(receipt_path.read_text())
if (len(receipt['checks'])!=21 or not all(receipt['checks'].values())
        or receipt['binary_sha256']!='feae0a48f737a7eab215f2e31be1d7938f45f894376dcdfda177f7e09f57d26b'
        or receipt['initial_sha256']!='8b561cdf05bfa06d2f42c5e7b782a2174719d5306967f4201b487d6113c4ab5b'
        or abs(receipt['time']-30.0)>1e-12):
    raise ValueError('Parent weighted t=30 physics gate failed')
source=Path(receipt['endpoint'])
if (source.parent!=parent or digest(source)!=receipt['endpoint_sha256']
        or abs(snapshot_time(source,55351)-30.0)>1e-10):
    raise ValueError('Parent weighted t=30 checkpoint changed or incomplete')
for label in ('control','candidate'):
    target=root+'/'+label+'/khout.000001.dat'
    if Path(target).exists():
        raise FileExistsError(target)
    shutil.copy2(source,target)
    if digest(Path(target))!=receipt['endpoint_sha256']:
        raise ValueError('Copied t=30 checkpoint differs from accepted parent')
with open(root+'/restart_origin.json','x') as out:
    json.dump({'receipt':str(receipt_path),'source':str(source),
               'sha256':receipt['endpoint_sha256'],'time':30.0},out,indent=2)
    out.write('\n')
PY_PARENT
export GFS_GPU_VOLW2=1 GFS_GPU_VOLW2_ITERS=64
export GFS_KEPLER_FIXED_POLYGON=1 GFS_KEPLER_HALO_CELLS=1
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0.06
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP=30.025 EUNHA_DUMP_DT=30.025 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS LAGMFM_GLASS_FILE EUNHA_KEPLER_INITIAL_AREAS_FILE
unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
unset GFS_GPU_WDOT_PROBE_MASK GFS_GPU_FACE_DUMP_CALL GFS_GPU_FACE_DUMP_IDS GFS_GPU_TESS_DUMP
nvidia-smi -L
for variant in control candidate; do
    cd "$root/$variant"
    sha256sum -c inputs.sha256
    export GFS_GPU_VOLW2_MASK="$PWD/disk_active_mask.bin"
    if [[ $variant == control ]]; then
        expected_binary=0e92efe71e2516b8a2f1b919ea6ba1c5c019532d1dd407331988b996dda46775
    else
        expected_binary=adeca5d738bce85a740b966fbbef04e3aaaba898ac7978eb5db797a620eaf5cd
    fi
    TIMEFORMAT='[RUN_TIME] wall=%R'
    { time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat 1 > log 2>&1; } 2> runtime.log
    /home/kjhan/miniconda3/bin/python3 audit_kepler_q2_weighted_equilibrium_20261027.py \
        initial.dat khout.000001.dat . --parent-time 30.0 --time 30.025 \
        --expected-binary-sha256 "$expected_binary" --output audit.json > audit.log 2>&1
done
cd "$root"
/home/kjhan/miniconda3/bin/python3 audit_kepler_weighted_fast_pair_20261027.py \
    control/audit.json candidate/audit.json \
    --control-binary-sha256 0e92efe71e2516b8a2f1b919ea6ba1c5c019532d1dd407331988b996dda46775 \
    --candidate-binary-sha256 adeca5d738bce85a740b966fbbef04e3aaaba898ac7978eb5db797a620eaf5cd \
    --parent-time 30.0 --time 30.025 --output pair_audit.json > pair_audit.log
