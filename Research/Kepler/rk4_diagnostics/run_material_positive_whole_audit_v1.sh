#!/usr/bin/env bash
# Run ONCE on the authorized CPU compute server after actual producer success.
# No polling, waiting daemon, rule proposal, GPU use or trajectory replay.
set -euo pipefail
root=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
producer=$root/material_full_positive_rule_v1_ada_mpi4
auditor=$root/material_positive_record_prefix_audit_v1_mpi4
out=${1:-$root/material_positive_rule_whole_audit_v1_mpi4}
test "$(hostname)" = LagEunha
test ! -e "$out"

# Check readiness BEFORE making any output directory. A live/partial producer
# is not a failed audit and is never adopted or restarted by this command.
python3 - "$root" "$producer" "$out" <<'PY'
import json, pathlib, sys
root,p,out=map(pathlib.Path,sys.argv[1:])
assert out.resolve().parent==root.resolve(), 'new audit output must be one child of this project'
assert not out.exists(), 'preserve previous audit/output'
for name in ['native_exit_code.txt','runner_exit_code.txt']:
    f=p/name
    assert f.is_file() and f.read_text().strip()=='0', 'WAIT: producer/runner not terminal successful'
c=json.loads((p/'contract.json').read_text())
assert c['candidate_ready_for_actual_RK_oracle'] and c['all_original_triangles_processed']
assert c['completed_triangles']==1138432 and c['original_positive_points']==2040256512
assert 0<c['retained_points_with_original_rule_fallbacks']<=2040256512//4
assert not c['production_rule_committed'] and not c['new_accepted_step']
assert abs(c['delta_anchor_U_over_U'])<1e-8 and abs(c['delta_anchor_Phi_over_U'])<1e-8
PY

mkdir "$out"
trap 'printf "WHOLE_RECORD_AUDIT_REJECTED_PRESERVE_OUTPUT\n" > "$out/status.txt"' ERR
source "$auditor/run_env.sh"
export I_MPI_CXX=g++ I_MPI_HYDRA_BOOTSTRAP=fork I_MPI_PIN=0 OMP_NUM_THREADS=1 OMP_PROC_BIND=false
unset OMP_PLACES GOMP_CPU_AFFINITY KMP_AFFINITY I_MPI_PMI_LIBRARY PMI_RANK PMI_SIZE PMI_FD PMI_PORT
sha256sum "$0" > "$out/runner.sha256"
sha256sum -c "$auditor/source.sha256" > "$out/auditor_source_integrity.txt"
(cd "$auditor" && sha256sum -c binary.sha256) > "$out/auditor_binary_integrity.txt"
sha256sum -c "$producer/source.sha256" > "$out/producer_source_integrity.txt"
sha256sum -c "$producer/binary.sha256" > "$out/producer_binary_integrity.txt"
sha256sum -c "$producer/input.sha256" > "$out/accepted_and_prepared_input_integrity.txt"
sha256sum -c "$producer/immutable_payload.sha256" > "$out/immutable_payload_integrity.txt"
sha256sum "$producer/contract.json" "$producer/native_exit_code.txt" "$producer/runner_exit_code.txt" \
    "$producer/run/positive_rule_r0.bin" "$producer/run/positive_rule_r1.bin" \
    "$producer/run/positive_rule_r2.bin" "$producer/run/positive_rule_r3.bin" > "$out/completed_rule_input.sha256"
mkdir "$out/run"
printf 'RUNNING_WHOLE_SAVED_RECORD_AUDIT_NO_PHYSICAL_EVOLUTION\n' > "$out/status.txt"
set +e
mpirun -np 4 "$auditor/audit.exe" "$root/material_full_initial_v2_h200_mpi4/run" "$producer" \
    "$root/material_full_cached_continuation_v1_ada_mpi4/run/continuous_step4/accepted_r" \
    "$root/material_full_resident_rhs_v1_ada_mpi4/run/next_RK_stage1/prepared_k1_r" \
    "$producer/immutable_payload.sha256" "$out/run" 0 > "$out/native.log" 2>&1
result=$?
set -e
printf '%s\n' "$result" > "$out/native_exit_code.txt"
test "$result" = 0
python3 - "$out" "$producer" <<'PY'
import json, pathlib, sys
d,p=map(pathlib.Path,sys.argv[1:]);r=json.loads((d/'run/contract.json').read_text());c=json.loads((p/'contract.json').read_text())
assert r['audit_pass'] and r['whole_rule_checked'] and not r['prefix_only']
assert r['triangles_checked']==1138432 and r['original_positive_points_accounted']==2040256512
assert r['selected_points_with_full_rule_fallbacks']==c['retained_points_with_original_rule_fallbacks']
assert not r['rule_adopted'] and not r['new_accepted_step'] and not r['long_Kepler_validated']
r.update(actual_native_exit=0, paired_RK4_oracle_still_required=True, no_field_queries_or_solve_replay=True)
(d/'contract.json').write_text(json.dumps(r,indent=2)+'\n')
PY
sha256sum -c "$out/completed_rule_input.sha256" > "$out/final_rule_integrity.txt"
printf 'WHOLE_RECORD_AUDIT_PASS_PAIRED_PHYSICS_STILL_REQUIRED\n' > "$out/status.txt"
