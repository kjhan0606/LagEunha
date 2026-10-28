#!/bin/bash
#SBATCH --job-name=kepler-inner-mass
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_InnerMass_%j.slurm.log

set -euo pipefail
module purge
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

SOURCE=/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs
INPUT=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026
READER=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
RUN=/gpfs/kjhan/LagForce/Kepler_InnerMass_Preflight_20261027
PYTHON=/home/kjhan/miniconda3/bin/python3

mkdir "$RUN"
cp "$SOURCE/audit_kepler_eulerian_inner_mass_20261027.py" "$RUN/"
cp "$SOURCE/plot_kepler_native_density_long_20261026.py" "$RUN/"

sha256sum \
  "$RUN/audit_kepler_eulerian_inner_mass_20261027.py" \
  "$RUN/plot_kepler_native_density_long_20261026.py" \
  "$READER/audit_snapshot.py" \
  "$READER/snapshot_offsets.json" \
  "$INPUT/q2/khout.000000.dat" \
  "$INPUT/q2/to0p825/khout.002700.dat" \
  "$INPUT/q2_profile_buffered_area_t001/khout.000000.dat" \
  "$INPUT/q2_profile_buffered_area_t001/khout.000001.dat" \
  "$INPUT/q2_profile_wide_area_t001/khout.000000.dat" \
  "$INPUT/q2_profile_wide_area_t001/khout.000001.dat" \
  > "$RUN/inputs.sha256"

"$PYTHON" "$RUN/audit_kepler_eulerian_inner_mass_20261027.py" \
  "$INPUT/q2/khout.000000.dat" \
  "$INPUT/q2/to0p825/khout.002700.dat" \
  --output "$RUN/full_q2_t020626.json" > "$RUN/full_q2.log" 2>&1

"$PYTHON" "$RUN/audit_kepler_eulerian_inner_mass_20261027.py" \
  "$INPUT/q2_profile_buffered_area_t001/khout.000000.dat" \
  "$INPUT/q2_profile_buffered_area_t001/khout.000001.dat" \
  --output "$RUN/sparse_narrow_t0001.json" > "$RUN/sparse_narrow.log" 2>&1

"$PYTHON" "$RUN/audit_kepler_eulerian_inner_mass_20261027.py" \
  "$INPUT/q2_profile_wide_area_t001/khout.000000.dat" \
  "$INPUT/q2_profile_wide_area_t001/khout.000001.dat" \
  --output "$RUN/sparse_wide_t0001.json" > "$RUN/sparse_wide.log" 2>&1

test -s "$RUN/full_q2_t020626.json"
test -s "$RUN/sparse_narrow_t0001.json"
test -s "$RUN/sparse_wide_t0001.json"
echo "Eulerian inner-mass preflight complete: $RUN"
