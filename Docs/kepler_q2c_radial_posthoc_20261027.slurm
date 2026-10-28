#!/usr/bin/env bash
#SBATCH --job-name=KeplerQ2cRadial
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
set -euo pipefail
module purge
run=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2c/halo12_after_failure_to0825
initial=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2c/khout.000000.dat
endpoint=$run/khout.006821.dat
output=$run/audit_radial_t0825.json
test -s "$initial"
test -s "$endpoint"
test ! -e "$output"
/home/kjhan/miniconda3/bin/python3 \
  /home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs/audit_kepler_native_radial_20261026.py \
  "$initial" "$endpoint" --mode baseline --output "$output" > "$run/audit_radial_posthoc.log"
sha256sum "$output"
