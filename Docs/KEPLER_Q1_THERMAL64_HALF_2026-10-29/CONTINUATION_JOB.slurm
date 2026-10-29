#!/bin/bash
#SBATCH --job-name=KepThermEarly0825
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --dependency=afterok:439903
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --chdir=/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64EarlyTo0825_H200_20261029
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64EarlyTo0825_H200_20261029/slurm-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_Q1_Thermal64EarlyTo0825_H200_20261029/slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=8 OMP_STACKSIZE=256M OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
hostname
nvidia-smi --query-gpu=name,uuid --format=csv,noheader
export PYTHONPATH=/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs
/home/kjhan/miniconda3/bin/python3 driver.py run
