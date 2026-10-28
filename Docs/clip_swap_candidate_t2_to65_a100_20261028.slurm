#!/usr/bin/env bash
# Isolated 384^2 W=0 candidate continuation from the exact q1 t=2 checkpoint.
# The production q1 trajectory remains untouched.

#SBATCH --job-name=KepClipSwapTo65
#SBATCH --partition=a100_pcie
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A100:1
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --chdir=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/diag_clip_swap_t2_to65_a100_20261028
#SBATCH --output=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/diag_clip_swap_t2_to65_a100_20261028/slurm-%j.out
#SBATCH --error=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/diag_clip_swap_t2_to65_a100_20261028/slurm-%j.err

set -euo pipefail

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2

project=/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026
q1="$project/q1"
run_dir="$q1/dteta_no_cap_from_t1p5_to65_20261028"
checkpoint="$run_dir/khout.016489.dat"
candidate_binary="$q1/diag_clip_swap_ab_t2_a100_pcie_20261028/candidate_binary.exe"
gate="$run_dir/audit_kepler_polar_long_gate_20261027.py"
radial="$run_dir/audit_kepler_native_radial_20261026.py"
python=/home/kjhan/miniconda3/bin/python3
mpi=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun
expected_candidate_sha=6718551061d70e24e09f0fe000a35b5cf1ac8c1a15b579f1e97e6d907b42baa0

for output in params.dat eunha2.exe khout.016489.dat khout.000000.dat log \
              run-manifest.txt audit_t6p5_gate.json audit_t6p5_radial.json; do
  if [[ -e "$output" || -L "$output" ]]; then
    echo "Refusing to overwrite $output" >&2
    exit 2
  fi
done
actual_candidate_sha=$(sha256sum "$candidate_binary" | awk '{print $1}')
[[ "$actual_candidate_sha" == "$expected_candidate_sha" ]] || {
  echo "candidate binary hash mismatch: $actual_candidate_sha" >&2; exit 2;
}

checkpoint_values=$("$python" - "$checkpoint" <<'PY'
import os
import struct
import sys

path = sys.argv[1]
size = os.path.getsize(path)
with open(path, "rb") as stream:
    magic, n, stride, stress = struct.unpack("<4i", stream.read(16))
    if (magic, n, stride, stress) != (-1, 147456, 576, 1):
        raise SystemExit("unexpected native checkpoint header")
    if size != 16 + n * stride + 52:
        raise SystemExit("incomplete native checkpoint")
    stream.seek(16 + n * stride)
    metadata = struct.unpack("<ddfiiddff", stream.read(52))
    if abs(metadata[0] - 2.0) > 1.0e-3:
        raise SystemExit(f"expected t=2 checkpoint, got t={metadata[0]}")
print(f"{metadata[0]:.15f} 16489")
PY
)
read -r start_time restart_step <<< "$checkpoint_values"

ln -s "$run_dir/params.dat" params.dat
ln -s "$candidate_binary" eunha2.exe
ln -s "$checkpoint" khout.016489.dat
ln -s "$q1/khout.000000.dat" khout.000000.dat

export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm
export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA=0
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP=6.5 EUNHA_DUMP_DT=0.5 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_KEPLER_HALO_CELLS=28
export GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=0.1
export GFS_GPU_LOCAL_PERIODIC_SYNC=0
unset GFS_GPU_TESS_PROFILE GFS_VARPRESS_GEOM_AUDIT GFS_GPU_SYNC_HALO_STATE
unset GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS
unset LAGMFM_GLASS_FILE EUNHA_KEPLER_INITIAL_AREAS_FILE

printf 'purpose=isolated candidate continuation\ncheckpoint=%s\ncheckpoint_sha256=%s\n' \
  "$checkpoint" "$(sha256sum "$checkpoint" | awk '{print $1}')" > run-manifest.txt
printf 'initial=%s\ninitial_sha256=%s\nstart_time=%s\nstop_time=6.5\nrestart_step=%s\n' \
  "$q1/khout.000000.dat" "$(sha256sum "$q1/khout.000000.dat" | awk '{print $1}')" \
  "$start_time" "$restart_step" >> run-manifest.txt
printf 'candidate_binary=%s\ncandidate_sha256=%s\nW=0\nmass_mode=fixed\n' \
  "$candidate_binary" "$actual_candidate_sha" >> run-manifest.txt

TIMEFORMAT='[RUN_TIME] wall=%R'
set +e
time "$mpi" -np 1 ./eunha2.exe params.dat "$restart_step" > log 2>&1
run_rc=$?
set -e
printf 'EXIT:%s\n' "$run_rc" >> log
printf 'exit=%s\n' "$run_rc" >> run-manifest.txt
[[ "$run_rc" -eq 0 ]] || exit "$run_rc"

"$python" "$gate" "$q1/khout.000000.dat" "$PWD" --time 6.5 \
  --output audit_t6p5_gate.json > audit_t6p5.log
endpoint=$("$python" -c 'import json; print(json.load(open("audit_t6p5_gate.json"))["endpoint"])')
"$python" "$radial" "$q1/khout.000000.dat" "$endpoint" --mode baseline \
  --output audit_t6p5_radial.json > audit_t6p5_radial.log
printf 'endpoint=%s\n' "$endpoint" >> run-manifest.txt
printf 'audit_t6p5_gate=PASS\naudit_t6p5_radial=PASS\n' >> run-manifest.txt
