#!/usr/bin/env bash
#SBATCH --job-name=KeplerA10Geom10
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
module purge
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [ "$#" -ne 2 ]; then
  echo 'Usage: runner RUN_DIRECTORY EXPECTED_COUNT' >&2
  exit 2
fi
/home/kjhan/miniconda3/bin/python3 - "$1" "$2" <<'PY'
import json
from pathlib import Path
import subprocess
import sys

run = Path(sys.argv[1]).resolve()
count = int(sys.argv[2])
if count not in (54312, 55931):
    raise ValueError("Unexpected sparse IC count")
physics = json.loads((run / "audit_t10.json").read_text())
state = json.loads((run / "audit_consistency_t10.json").read_text())
initial = Path(physics["initial"]).resolve()
endpoint = Path(physics["endpoint"]).resolve()
if (physics["n"] != count or state["n"] != count or
        abs(physics["time"] - 10.) > 1e-10 or
        abs(state["time"] - 10.) > 1e-10 or
        not all(physics["checks"].values()) or
        not all(state["checks"].values()) or
        initial != Path(state["initial"]).resolve() or
        endpoint != Path(state["endpoint"]).resolve() or
        physics["endpoint_sha256"] != state["endpoint_sha256"] or
        endpoint.parent != run):
    raise ValueError("The exact t=10 physical/state audit lineage has not passed")
output = run / "audit_t10_native_ballistic.json"
if output.exists():
    raise FileExistsError(output)
script = Path("/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs/") / \
    "audit_kepler_sparse_native_ballistic_20261027.py"
subprocess.run([sys.executable, str(script), str(initial), str(endpoint),
                "--expected-count", str(count), "--time", "10",
                "--output", str(output)], check=True)
PY
