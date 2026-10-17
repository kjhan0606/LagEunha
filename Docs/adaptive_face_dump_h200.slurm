#!/usr/bin/env bash
#SBATCH --job-name=LAG_FACE_ADAPT_256
#SBATCH --partition=h200
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:H200:2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/LagForce/adaptive_face_dump_%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/adaptive_face_dump_%j.err
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2

bundle=/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/adaptive_face_dump.bundle
baseline=/gpfs/kjhan/LagForce/Kepler_face_dump_Nx256_H200_halosync_call34188
runroot=/gpfs/kjhan/LagForce/Kepler_face_adaptive_256_${SLURM_JOB_ID}
test -s "$bundle"
test -s "$baseline/params.dat"
test ! -e "$runroot"
mkdir "$runroot"
git clone -q -b codex/adaptive-face-dump "$bundle" "$runroot/src"
git -C "$runroot/src" rev-parse HEAD > "$runroot/source_commit.txt"
printf '[BUILD] source=%s host=%s utc=%s\n' "$(cat "$runroot/source_commit.txt")" "$(hostname)" "$(date -u +%FT%TZ)"
make -C "$runroot/src" -j4 USE_CUDA=1 > "$runroot/build.log" 2>&1
test -x "$runroot/src/eunha2.exe"
cp "$runroot/src/eunha2.exe" "$runroot/eunha2.exe"
cp "$baseline/params.dat" "$runroot/params.dat"
sha256sum "$runroot/eunha2.exe" "$runroot/params.dat" > "$runroot/inputs.sha256"

export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
export SEDOV_PHASE1=1 EUNHA_IC=kepler
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export EUNHA_KEPLER_CX=12 EUNHA_KEPLER_CY=12
export EUNHA_KEPLER_RIN=2 EUNHA_KEPLER_ROUT=8
export EUNHA_KEPLER_EPS=0.047 EUNHA_KEPLER_DTETA=0.06
export HYDRO_TSTOP=7.5 EUNHA_DUMP_DT=1000
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1
export GFS_GPU_FACE_DUMP_HYD_THRESHOLD=1e-5
unset GFS_GPU_FACE_DUMP_CALL GFS_GPU_SYNC_HALO_STATE
unset GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_VOL_AUDIT GFS_FACE_AUDIT GFS_FACE_CHARGE_LIMIT GFS_GEOM_FACE_VEL
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT
unset SEDOV_LAGVOL SEDOV_WSMOOTH EUNHA_KEPLER_HR EUNHA_KEPLER_PROFILE

cd "$runroot"
printf '[PREFLIGHT] utc=%s cuda=%s\n' "$(date -u +%FT%TZ)" "${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi -L
set +e
/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 4 ./eunha2.exe params.dat > log 2>&1
rc=$?
set -e
printf 'EXIT:%s\n' "$rc" >> log
test "$rc" -eq 0
grep -q '\[FACE_DUMP_TRIGGER\]' log
shopt -s nullglob
rank0=(face_dump_r0_call*.tsv)
test "${#rank0[@]}" -eq 4
for file in "${rank0[@]}"; do
    call=${file#face_dump_r0_call}
    call=${call%.tsv}
    for rank in 0 1 2 3; do
        test -s "face_dump_r${rank}_call${call}.tsv"
    done
done
printf '[PREFLIGHT] complete utc=%s runroot=%s\n' "$(date -u +%FT%TZ)" "$runroot"
