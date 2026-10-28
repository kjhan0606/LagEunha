#!/usr/bin/env bash
#SBATCH --job-name=LagQ2Halo2755
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
set -euo pipefail
sha256sum -c inputs.sha256
python=/home/kjhan/miniconda3/bin/python3
initial=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_wide_area_t001/khout.000000.dat
checkpoint=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2_profile_wide_area_t001/to0p825/clipcompact_to6p5_a10/khdiag.step001522.start_geometry.dat
$python audit_kepler_periodic_halo_geometry_20261027.py "$initial" "$checkpoint" --image-margin 2.95 --halo-cells 32 --output halo32_t2755.json
