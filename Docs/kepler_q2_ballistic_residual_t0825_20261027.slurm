#!/bin/bash
#SBATCH --job-name=kepler-q2-residual
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_Q2_Residual_%j.slurm.log
set -euo pipefail
module purge
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

SOURCE=/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs
BASE=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2
NATIVE_RUN=$BASE/halo12_after_failure_to0825
INITIAL=$BASE/khout.000000.dat
Q1_BASE=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1
Q1_INITIAL=$Q1_BASE/khout.000000.dat
Q1_ENDPOINT=$Q1_BASE/halo12_from061876_to0825/khout.010856.dat
Q1_STRICT=$Q1_BASE/halo12_from061876_to0825/audit_long_t0825_strict.json
READER=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
RUN=/gpfs/kjhan/LagForce/Kepler_Q2_BallisticResidual_T0825_20261027
python=/home/kjhan/miniconda3/bin/python3

endpoint=$("$python" - "$SOURCE" "$NATIVE_RUN" <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from audit_kepler_polar_long_gate_20261027 import snapshot_time
paths = []
for path in Path(sys.argv[2]).glob("khout.*.dat"):
    if abs(snapshot_time(path) - .825) < 1e-10:
        paths.append(path)
if len(paths) != 1:
    raise SystemExit(f"Expected one complete q2 t=.825 native checkpoint; found {len(paths)}")
print(paths[0])
PY
)
test "$(sha256sum "$NATIVE_RUN/eunha2.exe" | cut -d' ' -f1)" = \
  7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284
test "$(sha256sum "$NATIVE_RUN/params.dat" | cut -d' ' -f1)" = \
  2898a04490002a24b766f0525426bfb41911e5733f54997b60673a5b5b69ef43
mkdir "$RUN"
cp "$SOURCE/audit_kepler_polar_long_gate_20261027.py" "$RUN/"
cp "$SOURCE/audit_kepler_polar_checkpoint_20261026.py" "$RUN/"
cp "$SOURCE/plot_kepler_polar_native_checkpoint_20261026.py" "$RUN/"
cp "$SOURCE/compare_kepler_polar_strict_t0825_20261027.py" "$RUN/"
sha256sum "$0" \
  "$RUN/audit_kepler_polar_long_gate_20261027.py" \
  "$RUN/audit_kepler_polar_checkpoint_20261026.py" \
  "$RUN/plot_kepler_polar_native_checkpoint_20261026.py" \
  "$RUN/compare_kepler_polar_strict_t0825_20261027.py" \
  "$READER/audit_snapshot.py" "$READER/snapshot_offsets.json" \
  "$NATIVE_RUN/eunha2.exe" "$NATIVE_RUN/params.dat" \
  "$INITIAL" "$endpoint" "$Q1_INITIAL" "$Q1_ENDPOINT" "$Q1_STRICT" > "$RUN/inputs.sha256"
sha256sum -c "$RUN/inputs.sha256"
"$python" "$RUN/audit_kepler_polar_long_gate_20261027.py" \
  "$INITIAL" "$NATIVE_RUN" --time .825 --output "$RUN/native_strict.json" \
  > "$RUN/native_strict.log" 2>&1
"$python" "$RUN/compare_kepler_polar_strict_t0825_20261027.py" \
  "$Q1_STRICT" "$RUN/native_strict.json" --output "$RUN/q1_q2_comparison.json" \
  > "$RUN/comparison.log" 2>&1
"$python" "$RUN/audit_kepler_polar_checkpoint_20261026.py" \
  "$INITIAL" "$endpoint" --output "$RUN/residual.json" \
  > "$RUN/residual.log" 2>&1
test -s "$RUN/residual.json"
"$python" "$RUN/plot_kepler_polar_native_checkpoint_20261026.py" \
  "$INITIAL" "$endpoint" --output "$RUN/q2_native_density.png" \
  > "$RUN/q2_plot.log" 2>&1
"$python" "$RUN/plot_kepler_polar_native_checkpoint_20261026.py" \
  "$Q1_INITIAL" "$Q1_ENDPOINT" --output "$RUN/q1_native_density.png" \
  > "$RUN/q1_plot.log" 2>&1
test -s "$RUN/q1_native_density.png"
test -s "$RUN/q2_native_density.png"
echo "q2 native/ballistic t=.825 decomposition complete: $RUN/residual.json"
