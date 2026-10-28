#!/usr/bin/env bash
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:A10:1
#SBATCH --mem=32G
#SBATCH --time=2:00:00
set -euo pipefail

run_dir=$(pwd)
python=/home/kjhan/miniconda3/bin/python3

for output in log audit.log audit_t1p1.json audit_radial_t1p1.json; do
	if [ -e "$output" ]; then
		echo "Refusing to overwrite $output" >&2
		exit 2
	fi
done

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0 cuda/13.0.2
export OMP_NUM_THREADS=2 OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}
export I_MPI_HYDRA_BOOTSTRAP=slurm

sha256sum -c inputs.sha256
restart_step=$($python -c 'import json; print(json.load(open("restart_plan.json"))["restart_step"])')
target_time=$($python -c 'import json; print(json.load(open("restart_plan.json"))["target_time"])')
eta=$($python -c 'import json; print(json.load(open("restart_plan.json"))["dt_eta"])')
initial=$($python -c 'import json; print(json.load(open("restart_plan.json"))["initial"])')

"$python" - <<'PY'
import json
import hashlib
import os
import struct

plan = json.load(open("restart_plan.json"))
path = plan["checkpoint"]

def digest(name):
	value = hashlib.sha256()
	with open(name, "rb") as stream:
		for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
			value.update(block)
	return value.hexdigest()

with open(path, "rb") as stream:
	header = struct.unpack("<4i", stream.read(16))
	stream.seek(-52, os.SEEK_END)
	time = struct.unpack("<d", stream.read(8))[0]
expected_size = 16 + 147456 * 576 + 52
assert header == (-1, 147456, 576, 1), header
assert os.path.getsize(path) == expected_size
assert abs(time - plan["restart_time"]) < 1e-12, (time, plan["restart_time"])
assert plan["target_time"] > time
assert plan["restart_step"] == 13200
for name, key in (("eunha2.exe", "binary_sha256"), ("params.dat", "params_sha256"),
			  (path, "checkpoint_sha256"), (plan["initial"], "initial_sha256")):
	assert digest(name) == plan[key], (name, key)
print(f"[RESTART_PREFLIGHT] time={time:.17g} step={plan['restart_step']} target={plan['target_time']:.17g} eta={plan['dt_eta']}")
PY

export SEDOV_PHASE1=1 EUNHA_IC=kepler EUNHA_KEPLER_PROFILE=hopkins
export EUNHA_KEPLER_CX=3 EUNHA_KEPLER_CY=3 EUNHA_KEPLER_RIN=0.5 EUNHA_KEPLER_ROUT=2
export EUNHA_KEPLER_FLOOR=0.01 EUNHA_KEPLER_WIDTH=0.1 EUNHA_KEPLER_EPS=0.01175 EUNHA_KEPLER_DTETA="$eta"
export EUNHA_DBG_LO=1000000000 EUNHA_DBG_HI=0
export HYDRO_TSTOP="$target_time" EUNHA_DUMP_DT=0.1 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
export GFS_RK4_FAILURE_SNAPSHOT=1 GFS_KEPLER_HALO_CELLS=28
export GFS_GPU_FINAL_VOLUME=1 GFS_GPU_RAW_FACES=1 GFS_VELOCITY_LIMITER=0
export GFS_FLOOR_LOG=1 GFS_E_AUDIT=1 GFS_VARPRESS_GEOM_CFL=.1 GFS_GPU_LOCAL_PERIODIC_SYNC=0
export GFS_GPU_TESS_PROFILE=1
unset GFS_DUMP_INITIAL EUNHA_KEPLER_MASS_FROM_CELL_VOLUME EUNHA_ALLOW_SPARSE_GLASS LAGMFM_GLASS_FILE EUNHA_KEPLER_INITIAL_AREAS_FILE
unset GFS_GPU_SYNC_HALO_STATE GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2
unset GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT GFS_RK4_TRANSACTION

TIMEFORMAT='[RUN_TIME] wall=%R'
set +e
time /opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin/mpirun -np 1 ./eunha2.exe params.dat "$restart_step" > log 2>&1
run_rc=$?
set -e
printf 'EXIT:%s\n' "$run_rc" >> log
if [ "$run_rc" -ne 0 ]; then
	exit "$run_rc"
fi

"$python" ./audit_kepler_polar_long_gate_20261027.py "$initial" "$run_dir" \
	--time "$target_time" --output audit_t1p1.json > audit.log
endpoint=$($python -c 'import json; print(json.load(open("audit_t1p1.json"))["endpoint"])')
"$python" ./audit_kepler_native_radial_20261026.py "$initial" "$endpoint" \
	--mode baseline --output audit_radial_t1p1.json > audit_radial.log
