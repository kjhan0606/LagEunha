#!/usr/bin/env bash
#SBATCH --job-name=LagQ2BWall825
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A10:1
#SBATCH --mem=32G
#SBATCH --time=00:20:00
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
{ hostname; nvidia-smi -L; printf 'CUDA_VISIBLE_DEVICES=%s\n' "${CUDA_VISIBLE_DEVICES:-unset}"; } > device.log
sha256sum -c inputs.sha256
unset LAGMFM_GLASS_FILE EUNHA_KEPLER_INITIAL_AREAS_FILE EUNHA_ALLOW_SPARSE_GLASS
unset EUNHA_KEPLER_MASS_FROM_CELL_VOLUME GFS_DUMP_INITIAL
export GFS_KEPLER_FIXED_POLYGON=1 GFS_KEPLER_HALO_CELLS=1
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0.06
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP=0.826 EUNHA_DUMP_DT=0.826 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
export GFS_GPU_TESS_PROFILE=1
export GFS_GPU_FACE_DUMP_CALL=1 GFS_GPU_FACE_DUMP_WALL_ONLY=1 GFS_GPU_TESS_DUMP=0
unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
TIMEFORMAT='[RUN_TIME] wall=%R'
set +e
time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat 1 > log 2>&1
run_rc=$?
set -e
printf 'EXIT:%s\n' "$run_rc" >> log
if [ "$run_rc" -ne 0 ]; then exit "$run_rc"; fi
set +e
/home/kjhan/miniconda3/bin/python3 audit_kepler_fixedpoly_wall_faces_20261027.py face_dump_r0_call1.tsv --allow-nonuniform --output audit_wall_call1.json > audit_wall.log 2>&1
wall_rc=$?
/home/kjhan/miniconda3/bin/python3 audit_kepler_bounded_short_replay_20261027.py . --time 0.826 --output audit_t0826.json > audit.log 2>&1
native_rc=$?
set -e
printf 'WALL_AUDIT_EXIT:%s NATIVE_AUDIT_EXIT:%s\n' "$wall_rc" "$native_rc" >> log
if [ "$wall_rc" -ne 0 ] || [ "$native_rc" -ne 0 ]; then exit 1; fi
