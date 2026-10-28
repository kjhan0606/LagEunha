#!/bin/bash
#SBATCH --job-name=kepler-q1-posthoc
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_Q1_Posthoc_%j.slurm.log

set -euo pipefail
module purge
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

SOURCE=/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs
BASE=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1
PARENT=$BASE/halo12_from061876_to0825
INITIAL=$BASE/khout.000000.dat
READER=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
RUN=/gpfs/kjhan/LagForce/Kepler_Q1_Halo12_Posthoc_20261027
PYTHON=/home/kjhan/miniconda3/bin/python3

test -s "$PARENT/audit_long_t0825.json"
ENDPOINT=$("$PYTHON" - "$PARENT/audit_long_t0825.json" "$PARENT" <<'PY'
import json, pathlib, sys
receipt=pathlib.Path(sys.argv[1])
parent=pathlib.Path(sys.argv[2]).resolve()
record=json.loads(receipt.read_text())
endpoint=pathlib.Path(record['endpoint']).resolve()
if abs(float(record['time'])-.825)>1e-10 or endpoint.parent!=parent:
    raise SystemExit('Unexpected q1 endpoint or time')
if not all(record['checks'].values()):
    raise SystemExit('Original q1 gate failed')
print(endpoint)
PY
)

mkdir "$RUN"
for name in \
  audit_kepler_polar_long_gate_20261027.py \
  audit_kepler_eulerian_inner_mass_20261027.py \
  audit_kepler_native_radial_20261026.py \
  audit_kepler_seam_checkpoint_20261026.py \
  plot_kepler_native_density_long_20261026.py; do
  cp "$SOURCE/$name" "$RUN/"
done
sha256sum \
  "$RUN"/*.py \
  "$READER/audit_snapshot.py" \
  "$READER/snapshot_offsets.json" \
  "$INITIAL" "$ENDPOINT" \
  "$PARENT/audit_long_t0825.json" \
  > "$RUN/inputs.sha256"

"$PYTHON" "$RUN/audit_kepler_polar_long_gate_20261027.py" \
  "$INITIAL" "$PARENT" --time .825 \
  --output "$RUN/strict_gate.json" > "$RUN/strict_gate.log" 2>&1
"$PYTHON" "$RUN/audit_kepler_native_radial_20261026.py" \
  "$INITIAL" "$ENDPOINT" --mode baseline \
  --output "$RUN/radial.json" > "$RUN/radial.log" 2>&1
"$PYTHON" "$RUN/audit_kepler_eulerian_inner_mass_20261027.py" \
  "$INITIAL" "$ENDPOINT" \
  --output "$RUN/inner_mass.json" > "$RUN/inner_mass.log" 2>&1

test -s "$RUN/strict_gate.json"
test -s "$RUN/radial.json"
test -s "$RUN/inner_mass.json"
echo "q1 posthoc audit complete: $RUN"
