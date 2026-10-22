#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
source "$root/run_env.sh"
export OMP_NUM_THREADS=1 I_MPI_HYDRA_BOOTSTRAP=fork
unset I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
cd "$root/$1"
test ! -e log
sha256sum ./eunha2.exe params.dat > inputs.sha256
python3 "$root/run_native.py"
