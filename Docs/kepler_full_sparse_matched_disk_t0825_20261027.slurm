#!/bin/bash
#SBATCH --job-name=kepler-full-sparse-match
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_FullSparseMatch_%j.slurm.log
set -euo pipefail
module purge
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
SOURCE=/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs
BASE=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026
RUN=/gpfs/kjhan/LagForce/Kepler_FullSparseMatchedDisk_T0825_20261027
FULL0=$BASE/q2/khout.000000.dat
FULL1=$BASE/q2/halo12_after_failure_to0825/khout.010856.dat
SPARSE0=$BASE/q2_profile_buffered_area_t001/khout.000000.dat
SPARSE1=$BASE/q2_profile_buffered_area_t001/to0p825/khout.000426.dat
READER=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
mkdir "$RUN"
cp "$SOURCE/audit_kepler_full_sparse_matched_disk_t0825_20261027.py" "$RUN/"
sha256sum "$0" "$RUN/audit_kepler_full_sparse_matched_disk_t0825_20261027.py" \
  "$READER/audit_snapshot.py" "$READER/snapshot_offsets.json" \
  "$FULL0" "$FULL1" "$SPARSE0" "$SPARSE1" > "$RUN/inputs.sha256"
sha256sum -c "$RUN/inputs.sha256"
/home/kjhan/miniconda3/bin/python3 \
  "$RUN/audit_kepler_full_sparse_matched_disk_t0825_20261027.py" \
  "$FULL0" "$FULL1" "$SPARSE0" "$SPARSE1" \
  --output "$RUN/result.json" > "$RUN/audit.log" 2>&1
test -s "$RUN/result.json"
echo "native full/sparse disk sites matched: $RUN/result.json"
