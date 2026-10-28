#!/usr/bin/env bash
#SBATCH --job-name=KepTess4GPU
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:H200:4
#SBATCH --mem=32G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_tess_4gpu_H200_20261026/slurm-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_tess_4gpu_H200_20261026/slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
cd /gpfs/kjhan/LagForce/Kepler_tess_4gpu_H200_20261026
sha256sum ./eunha2.exe > binary.sha256
nvidia-smi -L
export EUNHA_IC=kepler SEDOV_PHASE1=1
export EUNHA_KEPLER_CX=12 EUNHA_KEPLER_CY=12
export EUNHA_KEPLER_RIN=2 EUNHA_KEPLER_ROUT=8 EUNHA_KEPLER_EPS=0.047
export EUNHA_KEPLER_DTETA=0.06 HYDRO_TSTOP=0.003 EUNHA_DUMP_DT=1000
export GFS_GPU_TESS_PROFILE=1
set +e
/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 4 ./eunha2.exe params.dat > log 2>&1
rc=$?
set -e
printf 'EXIT:%s\n' "$rc" >> log
exit "$rc"
