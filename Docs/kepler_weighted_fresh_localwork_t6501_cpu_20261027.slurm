#!/usr/bin/env bash
#SBATCH --job-name=LagQ2LocalWork6501
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
early="$root/gpu_energy_neutral_capacity64_v5_fresh_localwork_t0826_ada_v1/audit_t0826.json"
late="$root/gpu_energy_neutral_capacity64_v5_fresh_localwork_t6501_ada_v1/audit_t6501.json"
"$python" audit_kepler_q2_weighted_fresh_localwork_t0826_20261027.py --receipt "$late" --start 6.5 --stop 6.501 --output local_work.json
"$python" audit_kepler_weighted_one_step_pdv_closure_20261027.py --early "$early" --late "$late" --output pdv_closure.json
