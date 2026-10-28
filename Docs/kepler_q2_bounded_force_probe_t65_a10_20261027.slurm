#!/usr/bin/env bash
#SBATCH --job-name=LagQ2Force65
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A10:1
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --nice=10000
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
export GFS_KEPLER_FIXED_POLYGON=1 GFS_KEPLER_HALO_CELLS=1
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0.06
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP=6.5005 EUNHA_DUMP_DT=6.5005 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
export GFS_GPU_FACE_DUMP_CALL=1 GFS_GPU_TESS_DUMP=0 GFS_GPU_FACE_DUMP_WALL_ONLY=0
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS LAGMFM_GLASS_FILE EUNHA_KEPLER_INITIAL_AREAS_FILE
unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
probe_root=$PWD
sha256sum -c "$probe_root/runner.sha256"
for variant in aligned jitter; do
    cd "$probe_root/$variant"
    sha256sum -c inputs.sha256
    read -r GFS_GPU_FACE_DUMP_IDS < selected_ids.txt
    export GFS_GPU_FACE_DUMP_IDS
    /home/kjhan/miniconda3/bin/python3 - <<'PY_RESTART'
from pathlib import Path
from prepare_kepler_sparse_area_restart_20261027 import snapshot_time
assert abs(snapshot_time(Path('khout.000001.dat'), 55351) - 6.5) < 1e-10
PY_RESTART
    TIMEFORMAT='[RUN_TIME] wall=%R'
    time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat 1 > log 2>&1
    test -s face_dump_r0_call1.tsv
    test -s face_input_r0_call1.tsv
    test -s face_input_meta_r0_call1.json
    /home/kjhan/miniconda3/bin/python3 audit_kepler_bounded_short_replay_20261027.py . \
        --time 6.5005 --output audit_short.json > audit_short.log
done
cd "$probe_root"
/home/kjhan/miniconda3/bin/python3 "$probe_root/audit_kepler_bounded_force_probe_t65_20261027.py" \
    "$probe_root" --output force_probe.json > force_probe.log
