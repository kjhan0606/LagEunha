#!/usr/bin/env bash
#SBATCH --job-name=KepCap16Audit
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --dependency=afterany:438420
#SBATCH --chdir=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/kernA_cap16_probe_t150_to1515_a100_20261028/extended_audit_20261028
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/kernA_cap16_probe_t150_to1515_a100_20261028/extended_audit_20261028/slurm-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/kernA_cap16_probe_t150_to1515_a100_20261028/extended_audit_20261028/slurm-%j.err
set -euo pipefail
module purge
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
sha256sum -c inputs.sha256
/home/kjhan/miniconda3/bin/python3 audit_kepler_cap16_pair_extended_20261028.py --output /home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs/KEPLER_CAP16_EXTENDED_AUDIT_2026-10-28.json
