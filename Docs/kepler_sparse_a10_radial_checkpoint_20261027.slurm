#!/usr/bin/env bash
#SBATCH --job-name=KeplerSparseRadial
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
set -euo pipefail
module purge
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [ "$#" -ne 4 ]; then
  echo 'usage: runner initial endpoint expected_time output' >&2
  exit 2
fi
initial=$1
endpoint=$2
expected_time=$3
output=$4
test -s "$initial"
test -s "$endpoint"
test ! -e "$output"
python=/home/kjhan/miniconda3/bin/python3
"$python" /home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs/audit_kepler_native_radial_20261026.py \
  "$initial" "$endpoint" --mode baseline --output "$output"
"$python" - "$output" "$expected_time" <<'PY'
import json
import math
import sys

receipt = json.load(open(sys.argv[1]))
expected = float(sys.argv[2])
if not math.isfinite(expected) or abs(receipt["time"] - expected) > 1e-10:
    raise SystemExit("Radial checkpoint time mismatch")
if not all(receipt["checks"].values()):
    raise SystemExit("Radial state gate failed")
print("AUDIT_OK", receipt["time"], receipt["cells"])
PY
sha256sum "$output"
