#!/bin/bash
#SBATCH --job-name=KepQ1RestartDt
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --dependency=afterok:439333
#SBATCH --kill-on-invalid-dep=yes
#SBATCH --chdir=/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/post_t10_restart_dt_h200_v1
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/post_t10_restart_dt_h200_v1/slurm-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/post_t10_restart_dt_h200_v1/slurm-%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=8 OMP_STACKSIZE=256M OPENBLAS_NUM_THREADS=1
hostname
nvidia-smi --query-gpu=name,uuid --format=csv,noheader
sha256sum -c bundle.sha256
/home/kjhan/miniconda3/bin/python3 /gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/post_t10_restart_dt_h200_v1/run_kepler_q1_restart_dt_bundle_20261029.py run /gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/post_t10_restart_dt_h200_v1
