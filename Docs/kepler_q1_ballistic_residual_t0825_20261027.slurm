#!/bin/bash
#SBATCH --job-name=kepler-q1-residual
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_Q1_Residual_%j.slurm.log

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
RUN=/gpfs/kjhan/LagForce/Kepler_Q1_BallisticResidual_T0825_20261027

mkdir "$RUN"
cp "$SOURCE/audit_kepler_polar_checkpoint_20261026.py" "$RUN/"
sha256sum \
  "$RUN/audit_kepler_polar_checkpoint_20261026.py" \
  "$READER/audit_snapshot.py" \
  "$READER/snapshot_offsets.json" \
  "$INITIAL" "$ENDPOINT" > "$RUN/inputs.sha256"

/home/kjhan/miniconda3/bin/python3 \
  "$RUN/audit_kepler_polar_checkpoint_20261026.py" \
  "$INITIAL" "$ENDPOINT" --output "$RUN/residual.json" \
  > "$RUN/audit.log" 2>&1
test -s "$RUN/residual.json"
echo "q1 native versus exact-orbit geometry audit complete: $RUN/residual.json"
