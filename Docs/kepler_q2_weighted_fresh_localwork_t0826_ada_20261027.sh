#!/usr/bin/env bash
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
unset I_MPI_HYDRA_BOOTSTRAP
sha256sum -c inputs.sha256
{ hostname; nvidia-smi -L; } > device.log
export GFS_GPU_VOLW2=1 GFS_GPU_VOLW2_ITERS=64
export GFS_GPU_VOLW2_MASK="$PWD/disk_active_mask.bin"
export GFS_KEPLER_FIXED_POLYGON=1 GFS_KEPLER_HALO_CELLS=1
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0.06
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP=0.826 EUNHA_DUMP_DT=0.001 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS LAGMFM_GLASS_FILE EUNHA_KEPLER_INITIAL_AREAS_FILE
unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
unset GFS_GPU_WDOT_PROBE_MASK GFS_GPU_FACE_DUMP_CALL GFS_GPU_FACE_DUMP_IDS GFS_GPU_TESS_DUMP
TIMEFORMAT='[RUN_TIME] wall=%R'
set +e
time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat 1 > log 2>&1
run_rc=$?
set -e
printf 'EXIT:%s\n' "$run_rc" >> log
if [ "$run_rc" -ne 0 ]; then exit "$run_rc"; fi
/home/kjhan/miniconda3/bin/python3 audit_kepler_q2_weighted_equilibrium_20261027.py \
  initial.dat khout.000001.dat . --parent-time 0.825 --time 0.826 \
  --expected-binary-sha256 feae0a48f737a7eab215f2e31be1d7938f45f894376dcdfda177f7e09f57d26b \
  --output audit_t0826.json > audit.log 2>&1
