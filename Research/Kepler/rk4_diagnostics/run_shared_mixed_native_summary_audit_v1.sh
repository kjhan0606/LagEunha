#!/usr/bin/env bash
if test "${GFS_CURVED_FROZEN_RUNNER:-0}" != 1; then exec bash /gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics/launch_curved_audit.sh "$0" "$@"; fi
set -euo pipefail
diag=/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/rk4_diagnostics
case "$(hostname)" in LagEunha|lageunha) ;; *) exit 2;; esac
cd "$diag"
dest=$diag/shared_mixed_native_summary_audit_v1
mkdir "$dest"
test "$(cat shared_mixed_native_points_audit_v2/native_exit_code.txt)" = 0
cp audit_shared_mixed_native_summary_v1.py "$dest/source.py"
sha256sum "$dest/source.py" audit_shared_atlas_witness_v1.py audit_shared_pending29_projection_v1.py audit_shared_regular_q_transfer_v1.py audit_shared_clipped_target_v1.py > "$dest/source.sha256"
printf '%s\n' "$GFS_CURVED_LAUNCH_ENVELOPE" > "$dest/runner_envelope.txt"
set +e
python3 - "$dest/source.py" "$diag/audit_shared_mixed_native_summary_v1.py" > "$dest/contract.json" 2> "$dest/stderr.log" <<'PY'
from pathlib import Path
import sys
frozen,location=sys.argv[1:]
sys.path.insert(0,str(Path(location).parent))
exec(compile(Path(frozen).read_bytes(),location,'exec'),{'__file__':location,'__name__':'__main__'})
PY
result=$?
set -e
printf '%s\n' "$result" > "$dest/native_exit_code.txt"
sha256sum -c "$dest/source.sha256" > "$dest/final_source_integrity.txt"
sha256sum "$dest/contract.json" > "$dest/output.sha256"
printf 'MIXED_SUMMARY_AUDIT_EXIT_%s\n' "$result" > "$dest/status.txt"
exit "$result"
