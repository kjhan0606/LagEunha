#!/usr/bin/env bash
# One-step, read-only Nsight Compute diagnostic for the late-state Voronoi
# kernel; also logs the RK4 geometric-CFL limiting cell. Submit only after
# the requested checkpoint exists in source_dir.
# Required sbatch export: PROFILE_CHECKPOINT_NAME.
# The production run and its files are never used as this job's working dir.

#SBATCH --job-name=KepKernelANCU3
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A10:1
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --chdir=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/local_refined_w0_r075_n420_20261027/diag_kernelA_ncu_t15_a10_mpiwrap_20261028
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/local_refined_w0_r075_n420_20261027/diag_kernelA_ncu_t15_a10_mpiwrap_20261028/slurm-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/local_refined_w0_r075_n420_20261027/diag_kernelA_ncu_t15_a10_mpiwrap_20261028/slurm-%j.err

set -euo pipefail

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2

source_dir=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/local_refined_w0_r075_n420_20261027/to2p221_a10
checkpoint_name=${PROFILE_CHECKPOINT_NAME:?set PROFILE_CHECKPOINT_NAME to a complete khout.NNNNNN.dat}

if [[ ! "$checkpoint_name" =~ ^khout\.[0-9]{6}\.dat$ ]]; then
  echo "Invalid checkpoint name: $checkpoint_name" >&2
  exit 2
fi
if [[ ! -s "$source_dir/$checkpoint_name" ]]; then
  echo "Checkpoint is missing or empty: $source_dir/$checkpoint_name" >&2
  exit 2
fi

snapshot_times=$(/home/kjhan/miniconda3/bin/python3 - "$source_dir/$checkpoint_name" <<'PY'
import os
import struct
import sys

path = sys.argv[1]
size = os.path.getsize(path)
with open(path, "rb") as stream:
    magic, n, stride, stress = struct.unpack("<4i", stream.read(16))
    if magic != -1 or stride != 576 or stress != 1 or size != 16 + n * stride + 52:
        raise SystemExit("incomplete or incompatible native checkpoint")
    stream.seek(16 + n * stride)
    meta = struct.unpack("<ddfiiddff", stream.read(52))
    if n != 176400 or abs(meta[0] - 1.5) > 1.0e-3:
        raise SystemExit(f"expected a 420^2 checkpoint near t=1.5; got n={n}, t={meta[0]}")
print(f"{meta[0]:.12f} {meta[0] + 1.0e-5:.12f}")
PY
)
read -r start_time stop_time <<< "$snapshot_times"

step_text=${checkpoint_name#khout.}
step_text=${step_text%.dat}
restart_step=$((10#$step_text))
printf -v restart_link 'khout.%06d.dat' "$restart_step"

for output in params.dat eunha2.exe khout.000000.dat "$restart_link" \
              kernelA_late.ncu-rep ncu-run.log profile-input.txt; do
  if [[ -e "$output" ]]; then
    echo "Refusing to overwrite $output" >&2
    exit 2
  fi
done

ln -s "$source_dir/params.dat" params.dat
ln -s "$source_dir/eunha2.exe" eunha2.exe
ln -s "$source_dir/khout.000000.dat" khout.000000.dat
ln -s "$source_dir/$checkpoint_name" "$restart_link"
printf 'source_checkpoint=%s start_time=%.12f stop_time=%s restart_step=%s\n' \
  "$checkpoint_name" "$start_time" "$stop_time" "$restart_step" > profile-input.txt

export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0.06
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP="$stop_time" EUNHA_DUMP_DT=0.5 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_KEPLER_HALO_CELLS=28
export GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
export GFS_GPU_TESS_PROFILE=0 GFS_VARPRESS_GEOM_AUDIT=1
unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS
unset LAGMFM_GLASS_FILE EUNHA_KEPLER_INITIAL_AREAS_FILE

set +e
/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 \
  /opt/ohpc/pub/cuda/13.0.2/bin/ncu \
  --config-file 0 --target-processes all \
  --kernel-name-base mangled --kernel-name 'regex:.*voronoi_tessellate_kernel.*' \
  --launch-count 1 --replay-mode kernel --clock-control none \
  --section LaunchStats --section Occupancy --section SpeedOfLight \
  --section WarpStateStats --section InstructionStats \
  --export kernelA_late \
  ./eunha2.exe params.dat "$restart_step" > ncu-run.log 2>&1
run_rc=$?
set -e
printf 'NCU_EXIT:%s\n' "$run_rc" >> ncu-run.log
exit "$run_rc"
