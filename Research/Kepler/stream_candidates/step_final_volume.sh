#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
source "$root/run_env.sh"
export OMP_NUM_THREADS=1 I_MPI_HYDRA_BOOTSTRAP=fork
export EUNHA_DUMP_DT=0.0001 GFS_GPU_FACE_DUMP_CALL=0
unset I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
for variant in "$@"; do
 case "$variant" in *gpu*) export GFS_GPU_FINAL_VOLUME=1;; *) export GFS_GPU_FINAL_VOLUME=0;; esac
 case "$variant" in *t8) export KEPLER_RESTART_STEP=13378 HYDRO_TSTOP=8.00014;; *) export KEPLER_RESTART_STEP=9570 HYDRO_TSTOP=6.0006;; esac
 cd "$root/$variant"
 test ! -e log
 sha256sum ./eunha2.exe params.dat > inputs.sha256
 python3 "$root/run_native.py"
done
