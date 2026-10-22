#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
source "$root/run_env.sh"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 "$root/rk4_diagnostics/run_ssprk3_cases.py" --phase fresh
