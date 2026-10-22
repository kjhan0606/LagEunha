#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
source "$root/run_env.sh"
export OMP_NUM_THREADS=1 I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
cd "$root/$1"
sha256sum ./eunha2.exe params.dat > inputs.sha256
python3 "$root/run_singleton.py"
