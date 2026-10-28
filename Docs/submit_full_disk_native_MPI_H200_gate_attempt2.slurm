#!/bin/bash
#SBATCH --job-name=LAGFORCE_FULLMOMENT_MPI_GATE
#SBATCH --partition=h200
#SBATCH --chdir=/gpfs/kjhan/LagForce/Codex_full_disk_native_MPI_H200_gate_20261019
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:2
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/LagForce/Codex_full_disk_native_MPI_H200_gate_20261019/slurm-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Codex_full_disk_native_MPI_H200_gate_20261019/slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb intel/compiler/2025.3.0 intel/mpi/2021.17 cuda/13.0.2
export OMP_NUM_THREADS=8 OMP_PROC_BIND=false OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 I_MPI_PIN_DOMAIN=omp
export LD_LIBRARY_PATH=/home/kjhan/local/lib:/opt/ohpc/pub/compiler/gcc/13.2.0/lib64:${LD_LIBRARY_PATH:-}
cd /gpfs/kjhan/LagForce/Codex_full_disk_native_MPI_H200_gate_20261019
python3 /gpfs/kjhan/LagForce/Codex_full_disk_native_MPI_H200_gate_20261019/run_gate.py
