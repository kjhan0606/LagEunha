#!/usr/bin/env bash
#SBATCH --job-name=LagQ2SnapshotCompare
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
root=/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027
w0="$root/jitter005_to65_a10_v1"
weighted="$root/gpu_energy_neutral_capacity64_v5_recover_to65_a10_v1"
initial="$root/jitter005_t001_a10_v1/khout.000000.dat"
step=002553
w0time=$($python -c 'import sys;from pathlib import Path;from prepare_kepler_sparse_area_restart_20261027 import snapshot_time;print(repr(snapshot_time(Path(sys.argv[1]),55351)))' "$w0/khout.$step.dat")
wtime=$($python -c 'import sys;from pathlib import Path;from prepare_kepler_sparse_area_restart_20261027 import snapshot_time;print(repr(snapshot_time(Path(sys.argv[1]),55351)))' "$weighted/khout.$step.dat")
"$python" audit_kepler_bounded_native_evolution_20261027.py "$initial" "$w0" --time "$w0time" --through-endpoint --expected-binary-sha256 3fc66a5fd36d41d6191a7c26c7186db18aefe6c990b826604fa7bea6f75cd55e --output w0_receipt.json
"$python" audit_kepler_q2_weighted_equilibrium_20261027.py "$weighted/initial.dat" "$weighted/khout.000001.dat" "$weighted" --time "$wtime" --parent-time .825 --expected-binary-sha256 feae0a48f737a7eab215f2e31be1d7938f45f894376dcdfda177f7e09f57d26b --output weighted_receipt.json
"$python" audit_kepler_bounded_entropy_sign_20261027.py w0_receipt.json --output w0_entropy.json
"$python" audit_kepler_weighted_entropy_sign_20261027.py weighted_receipt.json --output weighted_entropy.json
"$python" audit_kepler_weighted_vs_w0_fresh_20261027.py w0_receipt.json weighted_receipt.json w0_entropy.json weighted_entropy.json --output comparison.json
