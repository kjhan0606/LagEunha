#!/usr/bin/env bash
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
source "$root/run_env.sh"
export OMP_NUM_THREADS=1 I_MPI_HYDRA_BOOTSTRAP=fork NATIVE_MPI_RANKS=1
export EUNHA_IC=affine_probe HYDRO_TSTOP=0.000000000001 EUNHA_DUMP_DT=0.000000001 GFS_GPU_FACE_DUMP_CALL=1 KEPLER_RESTART_STEP=1 GFS_GPU_FINAL_VOLUME=1
unset I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
python3 "$root/prepare_irregular_affine.py"
for variant in rigid_regular_lsf rigid_irregular_lsf rigid_irregular_gg0 rigid_irregular_gg strain_irregular_lsf; do
 cd "$root/irregular_affine_probe/$variant"
 test ! -e log
 sha256sum ./eunha2.exe params.dat khout.000001.dat > inputs.sha256
 python3 "$root/run_native.py"
done
