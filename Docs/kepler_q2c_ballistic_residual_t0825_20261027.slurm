#!/usr/bin/env bash
#SBATCH --job-name=KeplerQ2cGeom
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2c/cellwise-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2c/cellwise-%j.err
set -euo pipefail
module purge
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

source_dir=/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs
base=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2c
initial=$base/khout.000000.dat
endpoint=$base/halo12_after_failure_to0825/khout.006821.dat
reader=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
run=/gpfs/kjhan/LagForce/Kepler_Q2c_CellwiseGeometry_T0825_20261027
python=/home/kjhan/miniconda3/bin/python3

test -s "$initial"
test -s "$endpoint"
mkdir "$run"
cp "$source_dir/audit_kepler_polar_checkpoint_20261026.py" "$run/"
sha256sum "$0" "$run/audit_kepler_polar_checkpoint_20261026.py" \
  "$reader/audit_snapshot.py" "$reader/snapshot_offsets.json" \
  "$initial" "$endpoint" > "$run/inputs.sha256"
sha256sum -c "$run/inputs.sha256"
"$python" "$run/audit_kepler_polar_checkpoint_20261026.py" \
  "$initial" "$endpoint" --output "$run/residual.json" > "$run/residual.log"
test -s "$run/residual.json"
sha256sum "$run/residual.json"
