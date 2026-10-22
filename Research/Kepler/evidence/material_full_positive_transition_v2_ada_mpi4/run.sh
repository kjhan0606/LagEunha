#!/usr/bin/env bash
set -euo pipefail
dest=$(cd -- "$(dirname -- "$0")" && pwd)
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
base=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020
producer=$root/material_full_positive_rule_v1_ada_mpi4
test "$(hostname)" = LagEunha
test ! -e "$dest/run"
source "$dest/run_env.sh"
export I_MPI_CXX=g++ I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 OMP_NUM_THREADS=8 OMP_PROC_BIND=false CUDA_VISIBLE_DEVICES=0
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
trap 'printf "PAIRED_TRANSITION_SETUP_OR_EXECUTION_REJECTED\n" > "$dest/status.txt"' ERR
sha256sum -c "$dest/source.sha256" > "$dest/prelaunch_source_integrity.txt"
sha256sum -c "$dest/runner_source.sha256" > "$dest/prelaunch_runner_integrity.txt"
sha256sum -c "$dest/binary.sha256" > "$dest/prelaunch_binary_integrity.txt"
sha256sum -c "$dest/retained_build_objects.sha256" > "$dest/retained_build_integrity.txt"
sha256sum -c "$dest/immutable_payload.sha256" > "$dest/immutable_payload_integrity.txt"
sha256sum -c "$dest/accepted_input.sha256" > "$dest/accepted_input_integrity.txt"
sha256sum -c "$dest/prepared_k1_input.sha256" > "$dest/prepared_k1_input_integrity.txt"
sha256sum -c "$producer/source.sha256" > "$dest/producer_source_integrity.txt"
sha256sum -c "$producer/input.sha256" > "$dest/producer_input_integrity.txt"
sha256sum -c "$producer/binary.sha256" > "$dest/producer_binary_integrity.txt"
python3 - "$producer" "$dest" <<'PY'
import hashlib,json,pathlib,sys
p,d=map(pathlib.Path,sys.argv[1:]); c=json.loads((p/'contract.json').read_text())
assert (p/'native_exit_code.txt').read_text().strip()=='0'
assert (p/'runner_exit_code.txt').read_text().strip()=='0'
assert c['candidate_ready_for_actual_RK_oracle'] is True and c['all_original_triangles_processed'] is True
assert c['completed_triangles']==1138432 and c['original_positive_points']==2040256512
assert 0<c['retained_points_with_original_rule_fallbacks']<=2040256512//4
assert not c['production_rule_committed'] and not c['new_accepted_step']
assert abs(c['delta_anchor_U_over_U'])<1e-8 and abs(c['delta_anchor_Phi_over_U'])<1e-8
files=[p/'contract.json',p/'native_exit_code.txt',p/'runner_exit_code.txt',p/'binary.sha256',p/'source.sha256']+[p/'run'/f'positive_rule_r{r}.bin' for r in range(4)]
with (d/'positive_rule_input.sha256').open('x') as out:
    for f in files:
        h=hashlib.sha256()
        with f.open('rb') as inp:
            for chunk in iter(lambda:inp.read(8*1024*1024),b''):h.update(chunk)
        out.write(h.hexdigest()+'  '+str(f)+'\n')
PY
active=$(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader)
if test -n "$active"; then printf '%s\n' "$active" > "$dest/device_busy.txt"; exit 3; fi
export GFS_MATERIAL_INITIAL_INPUT=$root/material_full_initial_v2_h200_mpi4/run
export GFS_MATERIAL_PREPARE_REFERENCE=1 GFS_MATERIAL_BOUNDARY=material_constant_Pext GFS_MATERIAL_EVOLVE=1
export GFS_MATERIAL_REFERENCE_INPUT=$root/material_full_reference_v2_h200_mpi4/run
export GFS_MATERIAL_IMMUTABLE_MANIFEST=$dest/immutable_payload.sha256
export GFS_MATERIAL_RESUME_PREFIX=$root/material_full_cached_continuation_v1_ada_mpi4/run/continuous_step4/accepted_r
export GFS_MATERIAL_DEVICE_MASS=1 GFS_MATERIAL_CONTINUOUS=1 GFS_MATERIAL_POINT_CACHE=1
export GFS_MATERIAL_QUADRATURE_ORACLE=1 GFS_POSITIVE_RULE_INPUT=$producer
export GFS_MATERIAL_PREPARED_NEXT_RHS=$root/material_full_resident_rhs_v1_ada_mpi4/run/next_RK_stage1
export GFS_MATERIAL_SAVE_SOURCE_CONTEXT=1 GFS_MATERIAL_TARGET_TIME=1 GFS_MATERIAL_WALL_SECONDS=86400
export GFS_GPU_RAW_FACES=1 GFS_GPU_FINAL_VOLUME=1 GFS_GPU_FACE_DUMP_CALL=0 GFS_GPU_TESS_DUMP=0
export HYDRO_TSTOP=10 EUNHA_DUMP_DT=100 GFS_EXACT_STOP=1 GFS_FAIL_NONPOSITIVE_U=1
unset GFS_MATERIAL_PREPARE_NEXT_RHS GFS_DG_PROJECT_ONLY GFS_MOMENT_INIT_ONLY GFS_DG_ENRICHED GFS_SSPRK3
mkdir "$dest/run"
cp "$base/core_refinement/fine_integrated/khout.000001.dat" "$dest/run/khout.000001.dat"
cp "$base/core_refinement/fine_integrated/params.dat" "$dest/run/params.dat"
df -B1 "$dest" > "$dest/free_space_at_launch.txt"
python3 - "$dest" <<'PY'
import json,pathlib,shutil,sys
d=pathlib.Path(sys.argv[1]); free=shutil.disk_usage(d).free
assert free>8*1024**3
(d/'output_policy.json').write_text(json.dumps({'run_purpose':'one actual paired RK4 transition; no long production release','checkpoints':'four accepted rank shards; on rejection four failed-stage and four rolled-back shards','source_context':'only per-source summaries at the actual low/high RHS stages, no point atlas','full_native_snapshots':0,'wrapper_wall_seconds':7200,'native_cost_projection_budget_seconds':86400,'planned_new_output_budget_bytes':4*1024**3,'free_bytes_at_launch':free,'previous_outputs_preserved':True},indent=2)+'\n')
PY
nvidia-smi --query-gpu=index,uuid,name,memory.total,memory.used --format=csv > "$dest/device.txt"
printf 'RUNNING_ONE_ACTUAL_PAIRED_RK4_TRANSITION\n' > "$dest/status.txt"
cd "$dest/run"
env | sort > "$dest/effective_environment.txt"
set +e
timeout --signal=TERM --kill-after=60s 7200s mpirun -np 4 "$dest/eunha2_positive_transition.exe" params.dat 1 > "$dest/native.log" 2>&1
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum -c "$dest/binary.sha256" > "$dest/final_binary_integrity.txt"
sha256sum -c "$dest/accepted_input.sha256" > "$dest/final_accepted_input_integrity.txt"
sha256sum -c "$dest/prepared_k1_input.sha256" > "$dest/final_prepared_k1_integrity.txt"
sha256sum -c "$dest/positive_rule_input.sha256" > "$dest/final_positive_rule_integrity.txt"
if test "$result" = 0; then
python3 - "$dest" <<'PY'
import json,pathlib,sys
d=pathlib.Path(sys.argv[1]); result=json.loads((d/'run/full_continuous_result.json').read_text())
assert result['new_accepted_steps']==1 and result['total_accepted_steps']==5
assert result['accepted_time']>result['time_before'] and result['no_replayed_accepted_interval']
s=d/'run/continuous_step5'; step=json.loads((s/'step.json').read_text());assert step['finite_step_gates_pass']
assert json.loads((s/'numerical_restart_roundtrip.json').read_text())['physical_and_numerical_ledgers_match']
for stage in range(1,5):assert json.loads((s/f'stage{stage}/quadrature_difference.json').read_text())['same_stage_quadrature_gate_pass']
assert json.loads((s/'endpoint/quadrature_difference.json').read_text())['endpoint_quadrature_gate_pass']
result.update(one_actual_paired_RK4_transition_pass=True,physical_long_term_heating_resolved=False,long_Kepler_validated=False)
(d/'contract.json').write_text(json.dumps(result,indent=2)+'\n')
PY
printf 'ACTUAL_PAIRED_STEP_PASS_LONG_PRODUCTION_STILL_UNPROVED\n' > "$dest/status.txt"
else
printf 'PAIRED_TRANSITION_REJECTED_NATIVE_EXIT_%s\n' "$result" > "$dest/status.txt"
fi
exit "$result"
