#!/bin/bash
#SBATCH --job-name=kepler-sparse-geometry
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_Sparse_LongGeometry_%j.slurm.log

set -euo pipefail
module purge
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

SOURCE=/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs
INPUT=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026
READER=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
RUN=/gpfs/kjhan/LagForce/Kepler_Sparse_LongGeometry_20261027

# mkdir refuses to reuse a previous audit or silently replace a result.
mkdir "$RUN"
cp "$SOURCE/audit_kepler_sparse_full_long_geometry_20261027.py" "$RUN/"
cp "$SOURCE/audit_kepler_polar_long_geometry_20261026.py" "$RUN/"
cp "$SOURCE/audit_kepler_polar_checkpoint_20261026.py" "$RUN/"
cp "$SOURCE/plot_kepler_native_density_long_20261026.py" "$RUN/"

sha256sum \
  "$RUN/audit_kepler_sparse_full_long_geometry_20261027.py" \
  "$RUN/audit_kepler_polar_long_geometry_20261026.py" \
  "$RUN/audit_kepler_polar_checkpoint_20261026.py" \
  "$RUN/plot_kepler_native_density_long_20261026.py" \
  "$READER/audit_snapshot.py" \
  "$READER/snapshot_offsets.json" \
  "$INPUT/q2/khout.000000.dat" \
  "$INPUT/q2_profile_buffered_area_t001/khout.000000.dat" \
  "$INPUT/q2_profile_wide_area_t001/khout.000000.dat" \
  > "$RUN/inputs.sha256"

"/home/kjhan/miniconda3/bin/python3" -c \
  'import sys, numpy, scipy; print(sys.version, numpy.__version__, scipy.__version__)' \
  > "$RUN/python_versions.txt"

"/home/kjhan/miniconda3/bin/python3" \
  "$RUN/audit_kepler_sparse_full_long_geometry_20261027.py" \
  "$INPUT/q2/khout.000000.dat" \
  "$INPUT/q2_profile_buffered_area_t001/khout.000000.dat" \
  "$INPUT/q2_profile_wide_area_t001/khout.000000.dat" \
  --output "$RUN/result.json" \
  > "$RUN/audit.log" 2>&1

test -s "$RUN/result.json"
echo "Geometry audit complete: $RUN/result.json"
