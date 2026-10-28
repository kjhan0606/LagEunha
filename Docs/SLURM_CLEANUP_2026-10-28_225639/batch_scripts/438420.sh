#!/usr/bin/env bash
#SBATCH --job-name=KepCap16Pair
#SBATCH --partition=a100_pcie
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A100:1
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --dependency=afterany:438172
#SBATCH --chdir=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/kernA_cap16_probe_t150_to1515_a100_20261028
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/kernA_cap16_probe_t150_to1515_a100_20261028/slurm-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/kernA_cap16_probe_t150_to1515_a100_20261028/slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP=1.515 EUNHA_DUMP_DT=0.5 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export EUNHA_GPU_DT_DIAG=1 EUNHA_GPU_FACE_AREA_AUDIT=1 GFS_GPU_TESS_PROFILE=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_KEPLER_HALO_CELLS=28
export GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS LAGMFM_GLASS_FILE
unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
unset GFS_VARPRESS_GEOM_AUDIT
sha256sum -c static_inputs.sha256
(cd source_overlay && sha256sum -c ../source_overlay.sha256)
if ldd ./eunha2.cap16.exe | grep 'not found'; then
  echo 'Unresolved shared library in diagnostic executable' >&2
  exit 3
fi
if [ -e log ] || [ -e area_cfl_audit.json ] || [ -e cap16_vs_cap48_a100_comparison.json ]; then
  echo 'Refusing to overwrite existing candidate audit output' >&2
  exit 2
fi
CONTROL_DIR=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/dt_areaavg_audit_t150_to1515_a100_20261028
if [ -e "$CONTROL_DIR/log" ] || [ -e "$CONTROL_DIR/area_cfl_audit.json" ]; then
  echo 'Refusing to overwrite existing 48-corner A100 control output' >&2
  exit 2
fi
sha256sum -c pair_audit.sha256
(cd "$CONTROL_DIR" && bash /gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/kernA_cap16_probe_t150_to1515_a100_20261028/run_area_audit_paired.sh)
TIMEFORMAT='[RUN_TIME] wall=%R'
SECONDS=0
set +e
time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.cap16.exe params.dat 13596 > log 2>&1
RC=$?
set -e
printf 'EXIT:%s\n[PAIR_RUN_TIME] wall_seconds=%s\n' "$RC" "$SECONDS" >> log
if [ "$RC" -ne 0 ]; then exit "$RC"; fi
/home/kjhan/miniconda3/bin/python3 ./audit_area_cfl.py
/home/kjhan/miniconda3/bin/python3 ./compare_cap16_reference.py
