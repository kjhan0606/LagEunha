#!/usr/bin/env bash
#SBATCH --job-name=KeplerMPI2
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:H200:2
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --chdir=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/mpi_scale_w0_t15_to155_h200_20261028/np2_after438363_20261028
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/mpi_scale_w0_t15_to155_h200_20261028/np2_after438363_20261028/slurm-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/mpi_scale_w0_t15_to155_h200_20261028/np2_after438363_20261028/slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
PY=/home/kjhan/miniconda3/bin/python3
sha256sum -c inputs.sha256
for f in log audit_t1p55_gate.json audit_t1p55_radial.json; do
  if [ -e "$f" ]; then echo "Refusing to overwrite $f" >&2; exit 2; fi
done
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP=1.55 EUNHA_DUMP_DT=0.05 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_KEPLER_HALO_CELLS=28 GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1
export GFS_VELOCITY_LIMITER=0 GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=.1
export GFS_GPU_LOCAL_PERIODIC_SYNC=0
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS LAGMFM_GLASS_FILE
unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION GFS_RK4_FAILURE_SNAPSHOT
TIMEFORMAT='[RUN_TIME] wall=%R'
set +e
time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 2 ./eunha2.exe params.dat 13596 > log 2>&1
RC=$?
set -e
printf 'EXIT:%s\n' "$RC" >> log
if [ "$RC" -ne 0 ]; then exit "$RC"; fi
"$PY" ./audit_kepler_polar_long_gate_20261027.py ./khout.000000.dat "$PWD" --time 1.55 --output audit_t1p55_gate.json > audit_t1p55.log
ENDPOINT=$("$PY" -c 'import json;print(json.load(open("audit_t1p55_gate.json"))["endpoint"])')
"$PY" ./audit_kepler_native_radial_20261026.py ./khout.000000.dat "$ENDPOINT" --mode baseline --output audit_t1p55_radial.json > audit_t1p55_radial.log
