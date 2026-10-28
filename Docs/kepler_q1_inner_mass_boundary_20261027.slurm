#!/bin/bash
#SBATCH --job-name=kepler-q1-inner
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_Q1_Inner_%j.slurm.log

set -euo pipefail
module purge
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

SOURCE=/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs
BASE=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1
INITIAL=$BASE/khout.000000.dat
ENDPOINT=$BASE/halo12_from061876_to0825/khout.010856.dat
READER=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
RUN=/gpfs/kjhan/LagForce/Kepler_Q1_InnerMass_Boundary_20261027

mkdir "$RUN"
cp "$SOURCE/audit_kepler_eulerian_inner_mass_20261027.py" "$RUN/"
cp "$SOURCE/plot_kepler_native_density_long_20261026.py" "$RUN/"
sha256sum \
  "$RUN/audit_kepler_eulerian_inner_mass_20261027.py" \
  "$RUN/plot_kepler_native_density_long_20261026.py" \
  "$READER/audit_snapshot.py" \
  "$READER/snapshot_offsets.json" \
  "$INITIAL" "$ENDPOINT" > "$RUN/inputs.sha256"

/home/kjhan/miniconda3/bin/python3 \
  "$RUN/audit_kepler_eulerian_inner_mass_20261027.py" \
  "$INITIAL" "$ENDPOINT" --output "$RUN/inner_mass.json" \
  > "$RUN/audit.log" 2>&1
test -s "$RUN/inner_mass.json"
echo "q1 boundary-aware inner mass audit complete: $RUN/inner_mass.json"
