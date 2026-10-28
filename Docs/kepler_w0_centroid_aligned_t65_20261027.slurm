#!/usr/bin/env bash
#SBATCH --job-name=LagQ2CentroidA65
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:0
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --nice=200
set -euo pipefail
module purge
module load gnu13/13.2.0
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
sha256sum -c inputs.sha256
/home/kjhan/miniconda3/bin/python3 audit_kepler_w0_centroid_kinematics_20261027.py aligned_receipt.json --output aligned_t65.json
