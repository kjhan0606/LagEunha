#!/bin/bash
#SBATCH --job-name=unified-gresho_128x128
#SBATCH --partition=h200
#SBATCH --gres=gpu:H200:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/rotating_shear_suite_v1/gresho_128x128
#SBATCH --output=/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/rotating_shear_suite_v1/gresho_128x128/slurm-%j.out
set -euo pipefail
python3 /gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/rotating_shear_suite_v1/run_case.py /gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/rotating_shear_suite_v1/gresho_128x128
