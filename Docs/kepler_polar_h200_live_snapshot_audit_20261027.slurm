#!/usr/bin/env bash
#SBATCH --job-name=LagQ2H200Snapshot
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:0
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --nice=200
set -euo pipefail
module purge
module load gnu13/13.2.0
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
sha256sum -c inputs.sha256
python=/home/kjhan/miniconda3/bin/python3
initial=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/khout.000000.dat
run=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825/to2p221_h28
snapshot="$run/khout.019926.dat"
time=$($python -c 'import sys;from pathlib import Path;from audit_kepler_polar_long_gate_20261027 import snapshot_time;print(repr(snapshot_time(Path(sys.argv[1]))))' "$snapshot")
sha256sum "$initial" "$snapshot" "$run/eunha2.exe" /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/audit_snapshot.py > data.sha256
"$python" audit_kepler_polar_long_gate_20261027.py "$initial" "$run" --time "$time" --output audit.json
"$python" audit_kepler_native_radial_20261026.py "$initial" "$snapshot" --mode baseline --output radial.json
sha256sum -c data.sha256
echo 'PASS: Snapshot gate and radial audit; cumulative ledger was read from the live run and may include later steps.'
