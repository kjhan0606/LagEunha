#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
source "$root/run_env.sh"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=0
active=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)
if [[ -n "$active" ]]; then
    printf 'Ada is busy; no native shadow run launched.\n' >&2
    exit 3
fi
python3 "$root/rk4_diagnostics/run_ale_shadow_v1.py" --phase "${1:?patches or real}"
