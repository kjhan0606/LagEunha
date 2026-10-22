#!/usr/bin/env bash
#SBATCH --job-name=KepRadialV1
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --array=0-1
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/radial-native-%A_%a.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/radial-native-%A_%a.err
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
envelope=${1:?frozen native launch envelope required}
sha256sum -c "$envelope/source.sha256"
source "$envelope/run_env.sh"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_PROC_BIND=false
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY
case "$SLURM_ARRAY_TASK_ID" in 0) label=uniform;; 1) label=graded;; *) exit 2;; esac
set +e
python3 "$envelope/run_radial_resolution_v1.py" --root "$root" --case "$label" --version v2
result=$?
set -e
printf '%s\n' "$result" > "$envelope/native_${label}_exit_code.txt"
exit "$result"
