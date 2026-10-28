#!/usr/bin/env bash
#SBATCH --job-name=KeplerQ2cSafety
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2c/safety-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2c/safety-%j.err
set -euo pipefail
module purge
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

source_dir=/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs
snapshot=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2c/khout.000000.dat
reader=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
run=/gpfs/kjhan/LagForce/Kepler_Q2c_PeriodicSafety_T0825_20261027
python=/home/kjhan/miniconda3/bin/python3

test -s "$snapshot"
mkdir "$run"
for name in audit_kepler_q2c_periodic_safety_20261027.py \
  audit_kepler_polar_long_geometry_20261026.py \
  audit_kepler_polar_checkpoint_20261026.py; do
  cp "$source_dir/$name" "$run/$name"
done
sha256sum "$0" "$run"/*.py "$reader/audit_snapshot.py" \
  "$reader/snapshot_offsets.json" "$snapshot" > "$run/inputs.sha256"
sha256sum -c "$run/inputs.sha256"
"$python" "$run/audit_kepler_q2c_periodic_safety_20261027.py" \
  "$snapshot" --output "$run/result.json" > "$run/audit.log"
test -s "$run/result.json"
sha256sum "$run/result.json"
