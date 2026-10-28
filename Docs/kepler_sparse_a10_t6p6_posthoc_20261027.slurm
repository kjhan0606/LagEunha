#!/usr/bin/env bash
#SBATCH --job-name=KeplerA10T66Audit
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
module purge
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [ "$#" -ne 3 ]; then
  echo 'Usage: runner RUN_DIRECTORY INITIAL_CHECKPOINT EXPECTED_COUNT' >&2
  exit 2
fi
/home/kjhan/miniconda3/bin/python3 - "$1" "$2" "$3" <<'PY'
import json
from pathlib import Path
import subprocess
import sys

run = Path(sys.argv[1]).resolve()
initial = Path(sys.argv[2]).resolve()
count = int(sys.argv[3])
if count not in (54312, 55931):
    raise ValueError("Unexpected sparse IC count")
scripts = Path("/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs")
sys.path.insert(0, str(scripts))
from prepare_kepler_sparse_area_restart_20261027 import snapshot_time

if abs(snapshot_time(initial, count)) > 1e-12:
    raise ValueError("The reference is not a complete native t=0 checkpoint")
candidates = []
for path in run.glob("khout.*.dat"):
    try:
        time = snapshot_time(path, count)
    except (OSError, ValueError):
        continue
    if abs(time - 6.6) < .02:
        candidates.append((path, time))
if len(candidates) != 1:
    raise ValueError(f"Expected one complete t≈6.6 checkpoint; found {candidates}")
endpoint, time = candidates[0]
radial = run / "audit_t6p6_radial.json"
geometry = run / "audit_t6p6_native_ballistic.json"
if radial.exists() or geometry.exists():
    raise FileExistsError("A t≈6.6 audit receipt already exists")
python = sys.executable
subprocess.run([python, str(scripts / "audit_kepler_native_radial_20261026.py"),
                str(initial), str(endpoint), "--mode", "baseline",
                "--output", str(radial)], check=True)
subprocess.run([python, str(scripts / "audit_kepler_sparse_native_ballistic_20261027.py"),
                str(initial), str(endpoint), "--expected-count", str(count),
                "--time", repr(time), "--output", str(geometry)], check=True)
r = json.loads(radial.read_text())
g = json.loads(geometry.read_text())
if (abs(r["time"] - time) > 1e-10 or abs(g["time"] - time) > 1e-10 or
        not all(r["checks"].values()) or not all(g["checks"].values())):
    raise ValueError("The t≈6.6 radial or geometry audit did not pass")
disk = g["bins"][-1]
print(json.dumps({"time": time, "n": count,
                  "native_density_abs": disk["native_density_change_abs_mass_weighted"],
                  "ballistic_density_abs": disk["ballistic_density_change_abs_mass_weighted"],
                  "native_ballistic_residual_abs": disk["native_minus_ballistic_density_abs_mass_weighted"]}))
PY
