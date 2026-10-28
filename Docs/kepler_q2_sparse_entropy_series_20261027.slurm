#!/usr/bin/env bash
#SBATCH --job-name=KeplerEntropySeries
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_Q2SparseEntropySeries_%j.log
set -euo pipefail
module purge
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
SOURCE=/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs
BASE=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026
RUN=/gpfs/kjhan/LagForce/Kepler_Q2SparseEntropySeries_20261027
READER=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
mkdir "$RUN"
cp "$SOURCE/audit_kepler_q2_sparse_entropy_series_20261027.py" \
   "$SOURCE/audit_kepler_sparse_area_t0825_20261027.py" \
   "$SOURCE/prepare_kepler_sparse_area_restart_20261027.py" "$RUN/"
full0="$BASE/q2/khout.000000.dat"
full1="$BASE/q2/to0p825/khout.002700.dat"
full2="$BASE/q2/to0p825/khout.005403.dat"
full3="$BASE/q2/to0p825/khout.008119.dat"
full4="$BASE/q2/halo12_after_failure_to0825/khout.010856.dat"
sparse0="$BASE/q2_profile_buffered_area_t001/khout.000000.dat"
sparse1="$BASE/q2_profile_buffered_area_t001/to0p825/khout.000107.dat"
sparse2="$BASE/q2_profile_buffered_area_t001/to0p825/khout.000214.dat"
sparse3="$BASE/q2_profile_buffered_area_t001/to0p825/khout.000320.dat"
sparse4="$BASE/q2_profile_buffered_area_t001/to0p825/khout.000426.dat"
sha256sum "$0" "$RUN"/*.py "$READER/audit_snapshot.py" \
  "$READER/snapshot_offsets.json" "$full0" "$full1" "$full2" "$full3" "$full4" \
  "$sparse0" "$sparse1" "$sparse2" "$sparse3" "$sparse4" > "$RUN/inputs.sha256"
sha256sum -c "$RUN/inputs.sha256"
/home/kjhan/miniconda3/bin/python3 \
  "$RUN/audit_kepler_q2_sparse_entropy_series_20261027.py" \
  --full-initial "$full0" --full "$full1" --full "$full2" --full "$full3" --full "$full4" \
  --sparse-initial "$sparse0" --sparse "$sparse1" --sparse "$sparse2" --sparse "$sparse3" --sparse "$sparse4" \
  --output "$RUN/result.json" > "$RUN/audit.log"
test -s "$RUN/result.json"
echo "completed native entropy series: $RUN/result.json"
