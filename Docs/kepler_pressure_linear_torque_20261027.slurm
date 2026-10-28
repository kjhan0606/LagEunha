#!/usr/bin/env bash
#SBATCH --job-name=LagQ2LinearPTorque
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
/home/kjhan/miniconda3/bin/python3 test_kepler_pressure_linear_torque_20261027.py
/home/kjhan/miniconda3/bin/python3 audit_kepler_pressure_linear_torque_20261027.py w0_receipt.json --output w0.json
/home/kjhan/miniconda3/bin/python3 audit_kepler_pressure_linear_torque_20261027.py weighted_receipt.json --output weighted.json
