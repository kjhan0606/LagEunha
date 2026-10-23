#!/usr/bin/env bash
#SBATCH --job-name=kepler-cell-capacity
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=1G
#SBATCH --time=00:02:00
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/material_cell_condensation_capacity_v1_h200/slurm-%j.out
set -euo pipefail
dest=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/material_cell_condensation_capacity_v1_h200
source "$dest/run_env.sh"
sha256sum -c "$dest/source.sha256" > "$dest/source_integrity.txt"
sha256sum -c "$dest/binary.sha256" > "$dest/binary_integrity.txt"
printf '%s\n' "$SLURM_JOB_ID" > "$dest/slurm_job_id.txt"
set +e
"$dest/probe" 106707056260 "$dest/capacity.json"
result=$?
set -e
printf '%s\n' "$result" > "$dest/exit_code.txt"
exit "$result"
