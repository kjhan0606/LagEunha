#!/bin/bash
#SBATCH --job-name=lag-kh256-t10
#SBATCH --partition=h100
#SBATCH --gres=gpu:H100:1
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --chdir=/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/kh_t10_approved_20261109_v1
#SBATCH --output=/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/kh_t10_approved_20261109_v1/n256/slurm_%j.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
srun --ntasks=1 --cpus-per-task=16 --cpu-bind=cores python3 /gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/kh_t10_approved_20261109_v1/launch.py 256
